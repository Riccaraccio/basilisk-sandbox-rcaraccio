/**
# Coupled diffusion of the gas species

This module solves the diffusion of all the gas species of one phase in one
multigrid system. The Fick correction is part of the operator, so the flux
of one species can depend implicitly on the gradients of the other species.

The equation for the species $i$ is
$$
\theta \frac{Y_i^{n+1} - Y_i^n}{\Delta t} = \nabla\cdot F_i + r_i
$$
with the flux
$$
F_i = G_i - Y_i U, \qquad U = \sum_k G_k
$$
`G_i` is the Fick flux of the species, with the sign of $\alpha\nabla a$ of
`poisson.h`. `U` is the corrective mass flux. With $\sum_i Y_i = 1$ the sum
of the fluxes $F_i$ is zero.

Two forms of `G_i` are available:

* Without `MWs`: $G_i = \alpha D_i \nabla Y_i$.
* With `MWs` (the form of `MOLAR_DIFFUSION`):
  $G_i = \alpha D_i (MW_i/MW) \nabla X_i$. The unknown stays $Y$. The
  solver replaces $X$ in each cell with its linearization about $Y^n$:
  $$
  \tilde{X}_i = \frac{MW}{MW_i}\left(Y_i - Y_i^n \sigma\right), \qquad
  \sigma = MW \sum_j \frac{Y_j - Y_j^n}{MW_j}
  $$
  `MW` is the mixture molecular weight of $Y^n$. When $Y = Y^n$, the flux is
  the X-gradient flux of `multicomponent-varprop.h`.

`alpha` is the face coefficient without the diffusivity, thus
$\rho_f f_s f_m$. `DList` holds the cell diffusivities. `theta` must include
`cm`, as for `diffusion()`.

## The linearization of the corrective flux

The product $Y_i U$ is linearized with both factors implicit:
$$
Y_i U \approx [Y_i^{n+1}]_f\, U^n + [Y_i^n]_f \left(U^{n+1} - U^n\right)
$$
The first term moves the species with the lagged corrective flux. It is the
term that limits the timestep in the explicit scheme. Here it is implicit, so
it has no Courant limit. The second term couples the species.

Caution: do not freeze the whole factor $Y_i$ at $Y^n$. Then the first term
becomes explicit, and the Courant limit of a trace species stays.

`[q]_f` is the face value. With `upwind` (the default), it is the value of the
upwind cell for the sign of $U^n$. Without `upwind`, it is the mean of the two
cells. The two choices keep $\sum_i [Y_i^n]_f = 1$. If $\sum Y^n = 1$, then
$\sum Y^{n+1} = 1$ is the solution of the summed system, to the solver
tolerance.

Caution: the central face value can make a species negative when the cell
Peclet number $|U| \Delta/(\alpha D_i)$ is more than 2. This can occur for a
heavy species in a light carrier. The system then loses its diagonal
dominance, and the relaxation diverges: in `test/coupled-diffusion.c` it
diverges at a corrective Courant number between 1.6 and 4.1. Use `upwind`
there.

## The cells without the phase

A cell can have a face with $f_s > 0$ and hold no phase: the phase-change
cases set $Y = 0$ where the phase fraction is below `F_ERR`. There $MW = 0$
and $\tilde{X} = 0$ for any $Y$, so the row of that cell has no diagonal
except $\theta/\Delta t$, and $\theta$ is the floor `F_ERR`. The neighbours
then push a flux into it, and its $Y$ grows without a bound. In
`run/fatehi-combustion.c` this made the solve of the pore gas diverge.

Thus a face carries a flux only when the two cells hold the phase,
$\sum Y^n > 0.5$. A cell without the phase keeps $Y^n$ (plus its source).
`cd_nmasked` gives the number of faces with $\alpha > 0$ that the last call
closed. The per-species `diffusion()` lets a small flux into such a cell.

## Memory

The solver uses $4N$ temporary cell fields for $N$ species: the frozen
$Y^n$, the right-hand side, and the correction and residual of `mg_solve()`.
It also writes the coarse levels of `DList` and `theta`. */

#include "poisson.h"

struct CoupledDiffusion {
  int n;
  scalar * Yn, * D;
  face vector alpha, un;
  scalar theta, MW;
  double * MWs;
  double dt;
  bool upwind;
};

int cd_nmasked = 0; // the faces that the last call closed, see above

static inline double cd_face (double qm, double qp, double u, bool upwind)
{
  if (upwind)
    return u > 0. ? qm : qp;
  return (qm + qp)/2.;
}

/**
`cd_sigma()` gives $\sigma$ of the linearization in each cell. It is zero
without `MWs`. */

static void cd_sigma (scalar * Yl, scalar sg, struct CoupledDiffusion * p)
{
  scalar MW = p->MW;
  foreach() {
    sg[] = 0.;
    if (p->MWs && MW[] > 0.) {
      int k = 0;
      for (scalar Y in Yl) {
        scalar Yn = p->Yn[k];
        sg[] += (Y[] - Yn[])/p->MWs[k];
        k++;
      }
      sg[] *= MW[];
    }
  }
}

/**
`cd_corrective()` gives the corrective flux $U = \sum_k G_k$ on the faces,
for the fields `Yl` with the $\sigma$ of `sg`. */

static void cd_corrective (scalar * Yl, scalar sg, face vector U,
                           struct CoupledDiffusion * p)
{
  face vector alpha = p->alpha;
  scalar MW = p->MW;
  foreach_face() {
    U.x[] = 0.;
    double MWf = face_value (MW, 0);
    int k = 0;
    for (scalar Y in Yl) {
      scalar Yn = p->Yn[k], Dk = p->D[k];
      double A = alpha.x[]*face_value (Dk, 0);
      if (p->MWs) {
        A = MWf > 0. ? A*p->MWs[k]/MWf : 0.;
        double Xm = MW[-1]*(Y[-1] - Yn[-1]*sg[-1])/p->MWs[k];
        double Xp = MW[]*(Y[] - Yn[]*sg[])/p->MWs[k];
        U.x[] += A*(Xp - Xm)/Delta;
      }
      else
        U.x[] += A*(Y[] - Y[-1])/Delta;
      k++;
    }
  }
}

/**
## The residual

The residual uses the full fields. The face flux of each species is
$F_i = G_i - [Y_i]_f U^n - [Y_i^n]_f (U - U^n)$. As in `poisson.h`, the
fluxes go into a face vector first, so the coarse and fine faces of a tree
stay conservative. */

static double cd_residual (scalar * al, scalar * bl, scalar * resl,
                           void * data)
{
  struct CoupledDiffusion * p = (struct CoupledDiffusion *) data;
  face vector alpha = p->alpha, un = p->un;
  scalar theta = p->theta, MW = p->MW;
  double maxres = 0.;

  scalar sg[];
  face vector U[], g[];
  cd_sigma (al, sg, p);
  cd_corrective (al, sg, U, p);

  for (int i = 0; i < p->n; i++) {
    scalar Y = al[i], b = bl[i], res = resl[i];
    scalar Yn = p->Yn[i], Di = p->D[i];
    foreach_face() {
      double A = alpha.x[]*face_value (Di, 0), G;
      if (p->MWs) {
        double MWf = face_value (MW, 0);
        A = MWf > 0. ? A*p->MWs[i]/MWf : 0.;
        double Xm = MW[-1]*(Y[-1] - Yn[-1]*sg[-1])/p->MWs[i];
        double Xp = MW[]*(Y[] - Yn[]*sg[])/p->MWs[i];
        G = A*(Xp - Xm)/Delta;
      }
      else
        G = A*(Y[] - Y[-1])/Delta;
      g.x[] = G
        - cd_face (Y[-1], Y[], un.x[], p->upwind)*un.x[]
        - cd_face (Yn[-1], Yn[], un.x[], p->upwind)*(U.x[] - un.x[]);
    }
    foreach (reduction(max:maxres), nowarning) {
      res[] = b[] + theta[]/p->dt*Y[];
      foreach_dimension()
        res[] -= (g.x[1] - g.x[])/Delta;
      if (fabs (res[]) > maxres)
        maxres = fabs (res[]);
    }
  }
  return maxres;
}

/**
## The relaxation

The relaxation is a block Gauss-Seidel: in each cell it updates all the
species together. It first computes the local defect $e = r - L(\delta Y)$
with the fluxes of the operator on level `l`. It then solves $M c = e$ and
adds $c$ to the correction. $M$ is an approximation of the diagonal block of
the operator: it takes the face values $[Y_i^n]_f$ at the cell value.

With $a'_i = \sum_f A_{i,f}/\Delta^2$ (times $MW/MW_i$ for the molar form),
$m_k = MW/MW_k$ and $y = Y^n$ of the cell, the block is
$$
M = -\mathrm{diag}(d) + (a' \circ y)\, m^T + y \left(a' - (a'\cdot y)\, m\right)^T
$$
with $d_i = a'_i + o + \theta/\Delta t$, where $o$ is the outflow part of
$U^n$. Without `MWs`, $m = 0$. $M$ is a diagonal plus a rank-2 term, so the
Woodbury identity solves it with a $2\times 2$ system. The cost is $O(N)$
for each cell, and no $N\times N$ matrix is stored. */

static void cd_relax (scalar * dal, scalar * resl, int l, void * data)
{
  struct CoupledDiffusion * p = (struct CoupledDiffusion *) data;
  face vector alpha = p->alpha, un = p->un;
  scalar theta = p->theta, MW = p->MW;
  int n = p->n;
  bool molar = (p->MWs != NULL), upwind = p->upwind;

  foreach_level_or_leaf (l, nowarning) {
    double e[n], ap[n], yn[n], mk[n];
    double Xm[n], X0[n], Xp[n], G0[n], G1[n], A0[n], A1[n];
    double idt = theta[]/p->dt, o = 0.;

    for (int i = 0; i < n; i++) {
      scalar da = dal[i], res = resl[i], Yn = p->Yn[i];
      e[i] = res[] + idt*da[];
      ap[i] = 0.;
      yn[i] = Yn[];
      mk[i] = (molar && MW[] > 0.) ? MW[]/p->MWs[i] : 0.;
    }

    foreach_dimension() {

      /**
      $\tilde{X}$ of the correction in the cell and in its two neighbours.
      The correction is homogeneous, so its $\sigma$ has no $Y^n$ part. */

      double sm = 0., s0 = 0., sp = 0.;
      if (molar) {
        for (int k = 0; k < n; k++) {
          scalar da = dal[k];
          sm += da[-1]/p->MWs[k];
          s0 += da[]/p->MWs[k];
          sp += da[1]/p->MWs[k];
        }
        sm *= MW[-1], s0 *= MW[], sp *= MW[1];
      }
      double MWf0 = face_value (MW, 0), MWf1 = face_value (MW, 1);
      double U0 = 0., U1 = 0.;
      for (int k = 0; k < n; k++) {
        scalar da = dal[k], Yn = p->Yn[k], Dk = p->D[k];
        A0[k] = alpha.x[]*face_value (Dk, 0);
        A1[k] = alpha.x[1]*face_value (Dk, 1);
        if (molar) {
          A0[k] = MWf0 > 0. ? A0[k]*p->MWs[k]/MWf0 : 0.;
          A1[k] = MWf1 > 0. ? A1[k]*p->MWs[k]/MWf1 : 0.;
          Xm[k] = MW[-1]*(da[-1] - Yn[-1]*sm)/p->MWs[k];
          X0[k] = MW[]*(da[] - Yn[]*s0)/p->MWs[k];
          Xp[k] = MW[1]*(da[1] - Yn[1]*sp)/p->MWs[k];
        }
        else
          Xm[k] = da[-1], X0[k] = da[], Xp[k] = da[1];
        G0[k] = A0[k]*(X0[k] - Xm[k])/Delta;
        G1[k] = A1[k]*(Xp[k] - X0[k])/Delta;
        U0 += G0[k], U1 += G1[k];
      }

      for (int i = 0; i < n; i++) {
        scalar da = dal[i], Yn = p->Yn[i];
        double F0 = G0[i]
          - cd_face (da[-1], da[], un.x[], upwind)*un.x[]
          - cd_face (Yn[-1], Yn[], un.x[], upwind)*U0;
        double F1 = G1[i]
          - cd_face (da[], da[1], un.x[1], upwind)*un.x[1]
          - cd_face (Yn[], Yn[1], un.x[1], upwind)*U1;
        e[i] -= (F1 - F0)/Delta;
        ap[i] += (A0[i] + A1[i])/sq(Delta);
      }
      o += (upwind ?
            max (un.x[1], 0.) + max (- un.x[], 0.) :
            (un.x[1] - un.x[])/2.)/Delta;
    }

    /**
    The Woodbury solve of $M c = e$ with $M = -\mathrm{diag}(d) + u_1 v_1^T
    + u_2 v_2^T$. */

    double ay = 0.;
    for (int i = 0; i < n; i++) {
      if (molar)
        ap[i] *= mk[i];
      ay += ap[i]*yn[i];
    }

    double K11 = 1., K12 = 0., K21 = 0., K22 = 1., r1 = 0., r2 = 0.;
    bool ok = true;
    for (int i = 0; i < n; i++) {
      double d = ap[i] + o + idt;
      if (d <= 0.) {
        ok = false;
        break;
      }
      double u1 = ap[i]*yn[i], u2 = yn[i];
      double v1 = mk[i], v2 = ap[i] - ay*mk[i];
      K11 -= v1*u1/d, K12 -= v1*u2/d;
      K21 -= v2*u1/d, K22 -= v2*u2/d;
      r1 -= v1*e[i]/d, r2 -= v2*e[i]/d;
    }
    double det = K11*K22 - K12*K21;
    if (ok && fabs (det) > 1e-12) {
      double z1 = (r1*K22 - K12*r2)/det, z2 = (K11*r2 - K21*r1)/det;
      for (int i = 0; i < n; i++) {
        scalar da = dal[i];
        double d = ap[i] + o + idt;
        da[] += (ap[i]*yn[i]*z1 + yn[i]*z2 - e[i])/d;
      }
    }
  }
}

/**
## User interface

The function solves the system and returns the statistics of `mg_solve()`.
`rList` is optional. When it is given, it must have one field for each
species. The function does not change `rList` or `theta`. */

trace
mgstats coupled_diffusion (scalar * YList, double dt, scalar * DList,
                           face vector alpha, scalar theta,
                           scalar * rList = NULL,
                           double * MWs = NULL,
                           bool upwind = true,
                           double tolerance = 0.,
                           int nrelax = 4,
                           int minlevel = 0)
{
  if (dt == 0.) {
    mgstats s = {0};
    return s;
  }

  int n = list_len (YList);
  assert (list_len (DList) == n);
  assert (!rList || list_len (rList) == n);

  scalar * Yn = list_clone (YList), * b = list_clone (YList);


  /**
  Caution: `list_clone()` copies the attributes of the species, and the
  species of the phase-change cases are VOF tracers. Two of those attributes
  broke the multigrid:

  * the stencil flag `s_restriction` ("the restriction is up-to-date"), so
    `restriction()` can skip the clone;
  * the restriction `restriction_volume_average`. With it, the coarse levels
    of `Yn`, and of the residual that `mg_solve()` clones from `b`, came out
    as NaN in `run/fatehi-combustion.c` (axisymmetric). The cause was not
    found.

  `b` gets the plain average, as `poisson.h` has for its residual. `Yn`
  holds the value of the phase, not the tracer form. If the species still
  have the functions of the VOF tracers, `Yn` gets the plain average and
  the bilinear prolongation. If the caller gave the species functions for
  phase values (`phase-prolongation.h`, as `multicomponent-varprop.h`
  does), `Yn` keeps them. `set_restriction()` also clears the stencil
  flags. */

#if TREE
  for (scalar s in Yn) {
    if (s.restriction == restriction_volume_average) {
      set_restriction (s, restriction_average);
      set_prolongation (s, refine_bilinear);
    }
    else
      set_restriction (s, s.restriction);
  }
  for (scalar s in b)
    set_restriction (s, restriction_average);
#endif
  scalar MW[], sg[], ph[];
  face vector un[], am[];

  foreach() {
    MW[] = 0.;
    ph[] = 0.;
    double s = 0.;
    for (int i = 0; i < n; i++) {
      scalar Y = YList[i], Y0 = Yn[i], bi = b[i];
      Y0[] = Y[];
      bi[] = - theta[]/dt*Y[];
      if (rList) {
        scalar r = rList[i];
        bi[] -= r[];
      }
      if (MWs)
        s += Y[]/MWs[i];
      ph[] += Y[];
    }
    if (MWs && s > 0.)
      MW[] = 1./s;
  }

  /**
  `am` is `alpha` with the faces closed next to a cell without the phase. */

  int nm = 0;
  foreach_face (reduction(+:nm)) {
    am.x[] = alpha.x[];
    if (alpha.x[] > 0. && (ph[] < 0.5 || ph[-1] < 0.5))
      am.x[] = 0., nm++;
  }
  cd_nmasked = nm;

  struct CoupledDiffusion p = {n, Yn, DList, am, un, theta, MW, MWs, dt,
                               upwind};

  /**
  $U^n$ is the corrective flux of $Y^n$, thus with $\sigma = 0$. */

  cd_sigma (Yn, sg, &p);
  cd_corrective (Yn, sg, un, &p);

  /**
  The relaxation reads these fields on every level. */

  restriction (Yn);
  restriction (DList);
  restriction ({MW, theta});
  restriction ((scalar *){am, un});

  mgstats s = mg_solve (YList, b, cd_residual, cd_relax, &p,
                        nrelax, NULL, max(1, minlevel),
                        tolerance ? tolerance : TOLERANCE);

  delete (Yn), free (Yn);
  delete (b), free (b);
  return s;
}
