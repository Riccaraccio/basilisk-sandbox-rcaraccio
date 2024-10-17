extern scalar omega;
extern scalar zeta;
extern double rhoS;

scalar psi[];
face vector ubf[];
mgstats mgpsf;

trace
mgstats project_sv (face vector ubf, scalar psi,
    double dt = 1.,
    int nrelax = 4)
{
  scalar prod[];
  foreach()
    prod[] = (omega[]*f[]*zeta[]/rhoS)/dt;

  mgstats mgp = poisson (psi, prod, tolerance = TOLERANCE/sq(dt),
      nrelax = nrelax);

  foreach_face()
    ubf.x[] = -dt*face_gradient_x (psi, 0);

  return mgp;
}

psi[right] = neumann (neumann_pressure(ghost));
psi[left]  = neumann (- neumann_pressure(0));

#if AXI
ubf.n[bottom] = 0.;
ubf.t[bottom] = dirichlet(0);
psi[top]    = neumann (neumann_pressure(ghost));
#else // !AXI
#  if dimension > 1
psi[top]    = neumann (neumann_pressure(ghost));
psi[bottom] = neumann (- neumann_pressure(0));
#  endif
#  if dimension > 2
psi[front]  = neumann (neumann_pressure(ghost));
psi[back]   = neumann (- neumann_pressure(0));
#  endif
#endif // !AXI
