/**
# Darcy and Forchheimer terms
This file implements the Darcy and Forchheimer terms for flow in porous media.
They allow to account for the resistance to flow due to the presence of a porous matrix.

In this step, we want to solve the following equation for the superficial
velocity $\mathbf{u}$:
$$
\frac{\partial \mathbf{u}}{\partial t} = - \frac{1}{\rho_g}
\left[ \frac{\mu_g\,\mathbf{u}}{\mathbf{K}} +
\rho_g\,F\,\frac{|\mathbf{u}|\,\mathbf{u}}{\sqrt{\mathbf{K}}} \right],
\qquad
F = \frac{1.75}{\sqrt{150\,\epsilon_g^3}} .
$$
$\mu_g$ is the viscosity of the gas, not an effective viscosity. With
$K = \epsilon^3 d^2/(150(1-\epsilon)^2)$ the steady state is the Ergun
equation, $-\nabla p = \mu\mathbf{u}/K + \rho F|\mathbf{u}|\mathbf{u}/\sqrt{K}$.

This is the formulation without `POROUS_ADVECTION`: the momentum equation has
the inertia $\rho$ and the pressure mobility $1/\rho$. The drag rate must carry
the same porosity factor as the mobility, so here it carries none.

Caution: `POROUS_ADVECTION` changes the mobility to $\epsilon/\rho$ in
`project_sf`, but not in `centered_gradient()`. With that flag the drag of this
file gives $K_{\rm eff} = K\epsilon$, and the cell velocity in the porous region
is $\mathbf{u}/\epsilon$. Do not use `POROUS_ADVECTION` with this file until the
mobility of that flag is moved into `alpha`.

Note that the permeability tensor **K** (Da in the code) can reach very small values.
Therefore, an implicit treatment of the Darcy and Forchheimer
terms must and is here implemented.

## The coupling with the pressure

Write $\lambda = f(A + B)$ for the drag rate of a cell and $c = \lambda\Delta t$.
Over one step, with the pressure gradient $G = \nabla p/\rho$ held constant,
the exact solution of $\partial_t \mathbf{u} = -\lambda\mathbf{u} - G$ is
$$
\mathbf{u}^{n+1} = e^{-c}\,\mathbf{u}^{*} - \Delta t\,\varphi(c)\,G,
\qquad
\varphi(c) = \frac{1 - e^{-c}}{c} .
$$
The `viscous_term` event applies the factor $e^{-c}$ to the predictor. The
`advection_term` event applies $\varphi(c)$ to the face mobility `alpha`, so
both projections of the step and `centered_gradient()` see the drag. In a
steady state this gives the Darcy law $\mathbf{u} = -G/\lambda$ for every
$\Delta t$.

Without the factor on `alpha`, the steady state is
$\mathbf{u} = -\Delta t\,G/(1 - e^{-c})$. For $c \gg 1$ that is
$-\Delta t\,G$: the pressure then depends on $\Delta t$ and does not depend on
**K**. Compile with `-DDARCY_PRESSURE_COUPLING=0` to get that old scheme back.

Caution: do not apply $e^{-c}$ to `alpha`. The steady state is then
$-\Delta t\,G\,e^{-c}/(1 - e^{-c})$, and for $c \gg 1$ the Poisson problem
becomes singular.

`alpha` stays a mobility, not a specific volume. No module reads it back as a
density. `gravity.h` reads it to turn a force into an acceleration, and that
acceleration then gets the factor $\varphi$, which the exact solution above
also gives to a constant force.

Extern variables defined elsewhere:

+ *porosity*: scalar field representing the porosity of the medium
+ *f*: scalar field representing the volume fraction of the presudo-phase
+ *rhoGv_S*, *muGv_S*: scalar fields for variable density and viscosity
+ *rhoG*, *muG*: double values for constant density and viscosity
*/

#ifndef F_ERR
  #define F_ERR 1e-10
#endif

extern scalar porosity;
extern scalar f;

#ifdef VARPROP
extern scalar rhoGv_S, muGv_S;
#else
extern double rhoG, muG;
#endif

/**
The  permeability tensor **Da** is defiened as a constant vector.
This is to allow for the possibility to have different values in each direction
to account for anisotropic porous media. Units: m^2
*/

coord Da = {1e-10, 1e-10};

#ifndef DARCY_PRESSURE_COUPLING
# define DARCY_PRESSURE_COUPLING 1
#endif

/**
## The drag rate

`darcy_lambda.x` holds $\lambda = f(A + B)$ for the direction `x`, with
$A = \mu/(K\rho)$ and $B = F|\mathbf{u}|/\sqrt{K}$. Both schemes below use
this function, so `DARCY_PRESSURE_COUPLING` changes only the coupling with the
pressure, not the drag.

`rhoGv_S` and `muGv_S` are the density and the viscosity of the pore gas.
Do not use the Brinkman viscosity $\mu/\epsilon$ here. A case can give that
value to the viscous term through `mu`, and it is a different quantity. */

vector darcy_lambda[];

static void darcy_cell_rate (void)
{
  foreach() {
    foreach_dimension()
      darcy_lambda.x[] = 0.;
    if (f[] > F_ERR) {
      double e = porosity[]/f[];
      double F = 1.75/sqrt (150.*cube (e));
      double Umag = norm(u);

      double muGh, rhoGh;
      #ifdef VARPROP
      muGh = muGv_S[];
      rhoGh = rhoGv_S[];
      #else
      muGh = muG;
      rhoGh = rhoG;
      #endif

      foreach_dimension() {
        double A = muGh/(Da.x*rhoGh);   // Darcy term
        double B = F*Umag/sqrt(Da.x);   // Forchheimer term
        darcy_lambda.x[] = (A + B)*f[];
      }
    }
  }
}

#if DARCY_PRESSURE_COUPLING

/**
The `advection_term` event fills `darcy_lambda` once per step. The
`viscous_term` event then reads the same values, so the predictor and the
mobility use the same $c$. The Forchheimer term uses the velocity at the
start of the step. */

/**
`alphad` is the mobility with the drag. `alpha_base` is the mobility that the
other modules set: `fm` by default, or `alphav` of `two-phase.h` and
`variable-properties.h`. Those modules write `alphav` again in each
`properties` event, and they never read `alpha` back, so `alpha` can point to
`alphad` for the whole run. */

face vector alphad[];
(const) face vector alpha_base = unityf;

event defaults (i = 0) {
  foreach_dimension() {
    darcy_lambda.x.nodump = true;
    alphad.x.nodump = true;
  }
}

/**
## The mobility event

Same-name events run in reverse order of the declaration, so this event runs
before the `advection_term` event of `centered.h` (and of `POROUS_ADVECTION`).
It runs after every `properties` event of the step. It therefore sets the
mobility for the projection of `pf`, the projection of `p`, and
`centered_gradient()`.

Caution: do not move this to a `properties` event. A `properties` event of
this file runs BEFORE the one of `variable-properties.h`. That event then
writes `alphav` again, and the change of the mobility has no effect.

The projection of `pf` in `advection_term` uses `dt/2`, but this factor uses
`dt`. That projection only predicts the face velocity for the advection. */

event advection_term (i++) {
  darcy_cell_rate();

  if (alpha.x.i != alphad.x.i) {
    alpha_base = alpha;
    alpha = alphad;
  }

  foreach_face() {
    double c = dt*face_value (darcy_lambda.x, 0);
    double phi = c > 0. ? -expm1 (-c)/c : 1.;
    alphad.x[] = alpha_base.x[]*phi;
  }
}

/**
## Viscous term event
This event runs BEFORE the viscous solve of `centered.h`, because same-name
events run in reverse order of the declaration. It damps the predictor with
the factor $e^{-c}$ of the rate from the `advection_term` event. */

event viscous_term (i++) {
  foreach()
    foreach_dimension()
      u.x[] *= exp(-darcy_lambda.x[]*dt);
}

#else // !DARCY_PRESSURE_COUPLING

/**
## Viscous term event (old scheme)
This event runs BEFORE the viscous solve of `centered.h`, because same-name
events run in reverse order of the declaration. It damps the velocity to
account for the Darcy and Forchheimer resistance. The projection does not see
the drag.
*/

event defaults (i = 0) {
  foreach_dimension()
    darcy_lambda.x.nodump = true;
}

event viscous_term (i++) {
  darcy_cell_rate();
  foreach()
    foreach_dimension()
      u.x[] *= exp(-darcy_lambda.x[]*dt);
}

#endif // DARCY_PRESSURE_COUPLING

/**
## Previous implementation (commented out)
Here we used the default acceleration field *a* to add the Darcy and Forchheimer
contributions explicitly. This appoach is more integrated with the existing framework but
less efficient due to the explicit nature of the terms.
*/

// Old approach, works but explicit
// event defaults (i = 0) {
//   if (is_constant(a.x)) {
//     a = new face vector;
//     foreach_face() {
//       a.x[] = 0.;
//       dimensional (a.x[] == Delta/sq(DT));
//     }
//   }
// }

// event acceleration (i++){
//   face vector av = a;
//   foreach_face() {
//     double ff = face_value(f, 0);
//     if (ff > F_ERR) {
//       double ef = face_value(porosity, 0);
//       double F  = 1.75/pow (150*pow (ef, 3), 0.5);
//       double Umag = norm (uf);

//       // Darcy contribution, weighted by the face fraction of the interface
//       av.x[] -= alpha.x[]/(fm.x[] + SEPS)* (muG*ef/Da.x) *uf.x[] *ff; 

//       // Forcheimer contribution
//       av.x[] -= alpha.x[]/(fm.x[] + SEPS)* (F*ef*rhoG/pow(Da.x,0.5)) *Umag  *uf.x[] *ff;
//     }
//   }
// }
