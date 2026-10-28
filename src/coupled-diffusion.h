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
int cd_nsplit = 0;  // the splits of a step since the start, see below

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
adds $c$ to the correction. $M$ is the exact diagonal block of the operator
of `cd_residual()`.

With $a_{i,f} = A_{i,f}/\Delta^2$ on the face $f$ (times $MW_i/MW_f$ for the
molar form), $a'_i = \sum_f a_{i,f}$, $m_k = MW/MW_k$ and $y = Y^n$ of the
cell, the block is
$$
M = -\mathrm{diag}(d) + (a' \circ m \circ y)\, m^T + \sum_f y_f\, w_f^T,
\qquad w_{f,j} = m_j \left(a_{j,f} - \sum_k a_{k,f}\, m_k y_k\right)
$$
with $d_i = a'_i m_i + o + \theta/\Delta t$, where $o$ is the outflow part
of $U^n$. $y_f$ is the face value $[Y^n]_f$ of the residual. Without `MWs`,
$m_k = 1$ in $d$ and $w_f$, and the second term is absent.

The faces where the cell is upwind have $y_f = y$, so they share one term.
Each other face adds one term with the $Y^n$ of the neighbour (half of it
for the central face value), unless that value is close to $y$ (see
`cd_add_face()`). $M$ is therefore a diagonal plus at most
$2 + 2d$ terms of rank 1 in $d$ dimensions, and the Woodbury identity solves
it with a small dense system. The cost is $O(N)$ for each cell, and no
$N\times N$ matrix is stored.

Caution: do not replace $y_f$ with the cell value $y$ in the block. The
first version did that. At a sharp step of a major species the two differ
by up to the height of the step. For a step of 0.98 at an angle to the grid
the relaxation then diverged for $Fo \geq 4$ and the solve stopped with a
SIGFPE (`test/coupled-sharp-front.c`). */

#define CD_MAXTERMS (2 + 2*dimension) // the cell, the molar term, the faces

/**
`cd_add_face()` adds the term of the face `A` (the $A_{k,f}$ of the face)
when the neighbour gives a part `1 - wown` of the face value: `wown` is 1
if the cell is upwind, 0 if the neighbour is, 1/2 for the central value.
`yb` is the $Y^n$ of the neighbour. `vsub` collects the parts that the
term of the cell must not count.

When the face value differs from the cell value by less than `CD_YFTOL` for
all the species, the face stays in the term of the cell, as in the first
version. The block is then not exact, but the smoother needs only an
approximation there, and a term less saves $O(m N)$ operations. With 88
species, the exact terms on all the faces cost 30 per cent more per cycle
and gave the same number of cycles on smooth fields. */

#define CD_YFTOL 1e-2

static void cd_add_face (int n, bool molar, double idx2, const double * A,
                         const double * mk, const double * yb,
                         const double * yn, double wown, double * vsub,
                         double (* u)[n], double (* v)[n], int * m)
{
  if (wown >= 1.)
    return;
  double dy = 0., sa = 0.;
  for (int j = 0; j < n; j++)
    dy = max (dy, fabs (yb[j] - yn[j])), sa += A[j];
  if ((1. - wown)*dy < CD_YFTOL || sa == 0.)
    return;
  double beta = 0.;
  if (molar)
    for (int k = 0; k < n; k++)
      beta += A[k]*idx2*mk[k]*yn[k];
  for (int j = 0; j < n; j++) {
    double w = molar ? mk[j]*(A[j]*idx2 - beta) : A[j]*idx2;
    u[*m][j] = yb[j], v[*m][j] = (1. - wown)*w;
    vsub[j] += (1. - wown)*w;
  }
  (*m)++;
}

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
    double yb[n], vsub[n], u[CD_MAXTERMS][n], v[CD_MAXTERMS][n];
    double idt = theta[]/p->dt, o = 0.;
    int m = molar ? 2 : 1; // 0: the faces where the cell is upwind, 1: molar

    for (int i = 0; i < n; i++) {
      scalar da = dal[i], res = resl[i], Yn = p->Yn[i];
      e[i] = res[] + idt*da[];
      ap[i] = 0.;
      yn[i] = Yn[];
      mk[i] = (molar && MW[] > 0.) ? MW[]/p->MWs[i] : 0.;
      vsub[i] = 0.;
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

      /**
      The terms $y_f w_f^T$ of the two faces. On the face 0 the cell is
      upwind when $U^n \leq 0$, on the face 1 when $U^n > 0$, as in
      `cd_face()`. */

      for (int side = 0; side <= 1; side++) {
        double uf = side ? un.x[1] : un.x[];
        double wown = !upwind ? 0.5 : side ? (uf > 0.) : (uf <= 0.);
        if (wown < 1.) {
          for (int j = 0; j < n; j++) {
            scalar Yn = p->Yn[j];
            yb[j] = side ? Yn[1] : Yn[-1];
          }
          cd_add_face (n, molar, 1./sq(Delta), side ? A1 : A0, mk, yb, yn,
                       wown, vsub, u, v, &m);
        }
      }
    }

    /**
    The Woodbury solve of $M c = e$ with $M = -\mathrm{diag}(d) +
    \sum_c u_c v_c^T$: $c = (U z - e)/d$, where $K z = r$ with
    $K = I - V^T D^{-1} U$ and $r = -V^T D^{-1} e$. */

    /**
    The term of the cell: the sum of $w_f$ over all the faces is
    $m \circ (a' - (a' \cdot (m \circ y)))$, from which the faces with a
    term of their own are taken away. */

    double ay = 0.;
    if (molar)
      for (int i = 0; i < n; i++)
        ay += ap[i]*mk[i]*yn[i];
    for (int i = 0; i < n; i++) {
      v[0][i] = (molar ? mk[i]*(ap[i] - ay) : ap[i]) - vsub[i];
      u[0][i] = yn[i];
      if (molar) {
        ap[i] *= mk[i];
        u[1][i] = ap[i]*yn[i], v[1][i] = mk[i];
      }
    }

    double K[m][m + 1], d[n], z[m];
    for (int a = 0; a < m; a++)
      for (int b = 0; b <= m; b++)
        K[a][b] = (a == b);
    bool ok = true;
    for (int i = 0; i < n && ok; i++) {
      d[i] = ap[i] + o + idt;
      if ((ok = (d[i] > 0.)))
        for (int a = 0; a < m; a++) {
          double va = v[a][i]/d[i];
          for (int b = 0; b < m; b++)
            K[a][b] -= va*u[b][i];
          K[a][m] -= va*e[i];
        }
    }

    /**
    Gauss elimination with partial pivoting. A pivot below `1e-12` of the
    largest entry of its column skips the update of the cell. */

    for (int c = 0; c < m && ok; c++) {
      int pv = c;
      double big = 0.;
      for (int a = c; a < m; a++) {
        if (fabs (K[a][c]) > fabs (K[pv][c]))
          pv = a;
        big = max (big, fabs (K[a][c]));
      }
      for (int a = 0; a < c; a++)
        big = max (big, fabs (K[a][c]));
      if (big == 0. || fabs (K[pv][c]) <= 1e-12*big)
        ok = false;
      else {
        if (pv != c)
          for (int b = 0; b <= m; b++)
            swap (double, K[c][b], K[pv][b]);
        for (int a = c + 1; a < m; a++) {
          double fct = K[a][c]/K[c][c];
          for (int b = c; b <= m; b++)
            K[a][b] -= fct*K[c][b];
        }
      }
    }

    if (ok) {
      for (int a = m - 1; a >= 0; a--) {
        double s = K[a][m];
        for (int b = a + 1; b < m; b++)
          s -= K[a][b]*z[b];
        z[a] = s/K[a][a];
      }
      for (int i = 0; i < n; i++) {
        scalar da = dal[i];
        double c = - e[i];
        for (int b = 0; b < m; b++)
          c += u[b][i]*z[b];
        da[] += c/d[i];
      }
    }
  }
}

/**
## The multigrid loop

`cd_mg()` is the loop of `mg_solve()`, with the same tuning of `nrelax`,
plus a test of divergence. It stops when the residual is not finite or
grows above `CD_DIVERGE` times its initial value. `*failed` is then true,
and it is also true when the loop ends above the tolerance. With `stop`
false, the loop runs as `mg_solve()` and never stops early. */

#define CD_DIVERGE 1e3

static mgstats cd_mg (scalar * a, scalar * b, struct CoupledDiffusion * p,
                      int nrelax, int minlevel, double tolerance,
                      bool stop, bool * failed)
{
  scalar * da = list_clone (a), * res = list_clone (b);
  for (int ib = 0; ib < nboundary; ib++)
    for (scalar s in da)
      s.boundary[ib] = s.boundary_homogeneous[ib];

  mgstats s = {0};
  s.nrelax = nrelax > 0 ? nrelax : 4;
  double resb = s.resb = s.resa = cd_residual (a, b, res, p);
  *failed = false;
  for (s.i = 0;
       s.i < NITERMAX && (s.i < NITERMIN || s.resa > tolerance);
       s.i++) {
    mg_cycle (a, res, da, cd_relax, p, s.nrelax, minlevel, grid->maxdepth);
    s.resa = cd_residual (a, b, res, p);
    if (stop && !(s.resa <= CD_DIVERGE*s.resb)) {
      *failed = true;
      break;
    }
    if (s.resa > tolerance) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
        s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
        s.nrelax--;
    }
    resb = s.resa;
  }
  s.minlevel = minlevel;
  if (s.resa > tolerance)
    *failed = true;
  delete (res), free (res);
  delete (da), free (da);
  return s;
}

/**
## One step

`cd_step()` solves the system for one step. With `stop`, a failed solve
gives back $Y^n$ and sets `*failed`. */

static mgstats cd_step (scalar * YList, double dt, scalar * DList,
                        face vector alpha, scalar theta, scalar * rList,
                        double * MWs, bool upwind, double tolerance,
                        int nrelax, int minlevel, bool stop, bool * failed)
{
  int n = list_len (YList);
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

  mgstats s = cd_mg (YList, b, &p, nrelax, max(1, minlevel), tolerance,
                     stop, failed);

  /**
  A failed solve gives back $Y^n$, so that the caller can split the step. */

  if (stop && *failed)
    foreach()
      for (int i = 0; i < n; i++) {
        scalar Y = YList[i], Y0 = Yn[i];
        Y[] = Y0[];
      }

  delete (Yn), free (Yn);
  delete (b), free (b);
  return s;
}

/**
## User interface

The function solves the system and returns the statistics of the solve.
`rList` is optional. When it is given, it must have one field for each
species. The function does not change `rList` or `theta`.

When a solve fails (see `cd_mg()`), the function does two steps of `dt/2`
instead, down to `dt/2^CD_MAXSPLIT`. At that depth it keeps the last
iterate and prints a warning, as `mg_solve()` does. The statistics then
add the cycles of all the solves, and `resa` is the one of the last solve.
The rate sources `rList` stay the same in each part of the step.

Caution: the split costs at least one failed solve, thus up to `NITERMAX`
cycles. It is a safety net for rare cells, not a substitute for a smaller
timestep. A solve fails, for example, at a discontinuous step of 0.98 of
a major species on a change of level, at $Fo \geq 16$
(`test/coupled-sharp-front.c`). */

#define CD_MAXSPLIT 4

static mgstats cd_split (scalar * YList, double dt, scalar * DList,
                         face vector alpha, scalar theta, scalar * rList,
                         double * MWs, bool upwind, double tolerance,
                         int nrelax, int minlevel, int depth)
{
  bool failed = false, last = (depth >= CD_MAXSPLIT);
  mgstats s = cd_step (YList, dt, DList, alpha, theta, rList, MWs, upwind,
                       tolerance, nrelax, minlevel, !last, &failed);
  if (!failed)
    return s;
  if (last) {
    scalar v = YList[0];
    fprintf (ferr, "coupled-diffusion.h: warning: convergence for %s not "
             "reached at dt %g after %d splits\n  res: %g nrelax: %d "
             "tolerance: %g\n", v.name, dt, depth, s.resa, s.nrelax,
             tolerance), fflush (ferr);
    return s;
  }
  cd_nsplit++;
  mgstats s1 = cd_split (YList, dt/2., DList, alpha, theta, rList, MWs,
                         upwind, tolerance, nrelax, minlevel, depth + 1);
  mgstats s2 = cd_split (YList, dt/2., DList, alpha, theta, rList, MWs,
                         upwind, tolerance, nrelax, minlevel, depth + 1);
  s2.i += s.i + s1.i;
  s2.resb = s.resb;
  s2.nrelax = max (s2.nrelax, max (s.nrelax, s1.nrelax));
  return s2;
}

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
  assert (list_len (DList) == list_len (YList));
  assert (!rList || list_len (rList) == list_len (YList));
  return cd_split (YList, dt, DList, alpha, theta, rList, MWs, upwind,
                   tolerance ? tolerance : TOLERANCE, nrelax, minlevel, 0);
}
