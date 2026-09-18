/**
# Budget of the divergence after the projection

`project_sf()` in `navier-stokes/centered-phasechange.h` solves

$$
\sum_d \frac{u_{f,d}^{+} - u_{f,d}^{-}}{\Delta} = c_m\,\nabla\cdot\mathbf{u}
 = -\,\text{div\_source} ,
$$

with `div_source` = `gas_source` + `drhodt` (filtered with
`GAS_SOURCE_EXACT`). This header measures how well the face velocity `uf`
satisfies that identity.

A check that reads `gas_source` alone does not measure the solver. Without
`NO_EXPANSION` its residual is approximately $|\text{drhodt}|$. Read
`div_source`, because it is the source that the projection used.

The event below runs at `end_timestep`, after the `projection` event and
before `adapt`. After `adapt`, the cells that the grid refined or coarsened
do not satisfy the identity exactly. Thus a probe at the start of the next
step does not measure the solver either.

The event fills these values at each step:

* `divb_Qds`: the integral of `div_source` over the domain,
* `divb_Qdiv`: the integral of the discrete divergence $c_m\nabla\cdot\mathbf{u}$,
* `divb_resmax`: the maximum of $|c_m\nabla\cdot\mathbf{u} + \text{div\_source}|$,
* `divb_dsmax`: the maximum of $|\text{div\_source}|$, the scale of the residual.

The integrals use `sq(Delta)`, because the discrete divergence already
carries `cm[]`. After the projection, `divb_Qdiv = -divb_Qds` to the
tolerance of the solver. The multigrid residual `mgp.resa` is the residual of
the same equation divided by `dt`, so `divb_resmax` is approximately
`dt*mgp.resa`.

Include this header after `navier-stokes/centered-phasechange.h`. */

double divb_Qds = 0., divb_Qdiv = 0., divb_resmax = 0., divb_dsmax = 0.;

event end_timestep (i++) {
  double Qds = 0., Qdiv = 0., resmax = 0., dsmax = 0.;
  foreach (reduction(+:Qds) reduction(+:Qdiv)
           reduction(max:resmax) reduction(max:dsmax)) {
    double d = 0.;
    foreach_dimension()
      d += uf.x[1] - uf.x[];
    d /= Delta;
    Qds += div_source[]*sq(Delta);
    Qdiv += d*sq(Delta);
    resmax = max (resmax, fabs (d + div_source[]));
    dsmax = max (dsmax, fabs (div_source[]));
  }
  divb_Qds = Qds, divb_Qdiv = Qdiv, divb_resmax = resmax, divb_dsmax = dsmax;
}
