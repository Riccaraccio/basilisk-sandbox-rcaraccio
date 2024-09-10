scalar ps[];
face vector ufs[], ufext[];
mgstats mgpsf;

scalar zeta[];
scalar omega[];
trace
mgstats project_sv (face vector uf, scalar p,
     (const) face vector alpha = unityf,
     double dt = 1.,
     int nrelax = 4)
{
  scalar div[];
  foreach()
    div[] = (omega[]*f[])/dt;

  mgstats mgp = poisson (ps, div, alpha,
      tolerance = TOLERANCE/sq(dt), nrelax = nrelax);

  foreach_face()
    ufs.x[] = -dt*alpha.x[]*face_gradient_x (ps, 0);

  return mgp;
}

ps[right] = neumann (neumann_pressure(ghost));
ps[left]  = neumann (- neumann_pressure(0));

#if AXI
ufs.n[bottom] = 0.;
ufs.t[bottom] = dirichlet(0);
ps[top]    = neumann (neumann_pressure(ghost));
#else // !AXI
#  if dimension > 1
ps[top]    = neumann (neumann_pressure(ghost));
ps[bottom] = neumann (- neumann_pressure(0));
#  endif
#  if dimension > 2
ps[front]  = neumann (neumann_pressure(ghost));
ps[back]   = neumann (- neumann_pressure(0));
#  endif
#endif // !AXI

event end_timestep (i++, last)
{

  mgpsf = project_sv (ufs, ps, alpha, dt, mgpsf.nrelax);

}


