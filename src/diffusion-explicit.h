trace
int diffusion_explicit (scalar f, double dt,
    face vector D,
    (const) scalar r = zeroc,     // default 0
    (const) scalar beta = zeroc,  // default 0
    (const) scalar theta = unity) // default 1
{

  /**
  If *dt* is zero we don't do anything. */

  if (dt == 0.)
    return 1;

  /**
  Otherwise we solve the diffusion equation. */

  face vector J[];
  foreach_face()
    J.x[] = D.x[]*face_gradient_x (f, 0);

  foreach() {
    f[] += dt/theta[]*beta[]*f[];       // add "implicit" source
    f[] += dt/theta[]*r[];              // add explicit source
    foreach_dimension()
      f[] += dt*(J.x[1] - J.x[])/(Delta*theta[]);
  }

  return 0;
}

trace
int diffusion_explicit_molar (scalar f, scalar fx, double dt,
    face vector D,
    (const) scalar r = zeroc,     // default 0
    (const) scalar beta = zeroc,  // default 0
    (const) scalar theta = unity) // default 1
{

  /**
  If *dt* is zero we don't do anything. */

  if (dt == 0.)
    return 1;

  /**
  Otherwise we solve the diffusion equation. */

  face vector J[];
  foreach_face()
    J.x[] = D.x[]*face_gradient_x (fx, 0);

  foreach() {
    f[] += dt/theta[]*beta[]*f[];       // add "implicit" source
    f[] += dt/theta[]*r[];              // add explicit source
    foreach_dimension()
      f[] += dt*(J.x[1] - J.x[])/(Delta*theta[]);
  }

  return 0;
}

