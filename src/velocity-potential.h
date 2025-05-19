extern scalar omega;
extern scalar zeta;
extern double rhoS;

#include "diffusion.h"

bool shift_prod = false;
bool diffuse_prod = false;

scalar psi[];
face vector ubf[];
mgstats mgpsf;
double TOLERANCE_SOLID = 1e-5;
extern face vector fsS;

void shift_field (scalar fts, scalar f, int dir);

scalar fts[];
scalar prod[];
extern scalar fS;

trace
mgstats project_sv (face vector ubf, scalar psi,
    (const) face vector alpha = unityf,
    int nrelax = 4)
{
  foreach()
    prod[] = omega[]*f[]*zeta[]*cm[]/rhoS;
  
  // double intbefore = 0.;
  // foreach()
  //   intbefore += prod[];

  if (shift_prod) 
    shift_field (prod, f, 1);
  
  if (diffuse_prod) {
    face vector D[];
    scalar theta[];

    foreach()
      theta[] = cm[]*max(fS[], F_ERR);

    foreach_face()
      D.x[] = fm.x[]*fsS.x[]*1e-5;

    diffusion (prod, dt, D=D, theta=theta);
  }

  // double intafter = 0.;
  // foreach()
  //   intafter += prod[];

  // double intdiff = intbefore > F_ERR ? (intbefore-intafter) /intbefore : 0.;
  // fprintf (stderr, "intdiff = %g\n", intdiff);

  mgstats mgp = poisson (psi, prod, alpha,
      tolerance = TOLERANCE_SOLID, nrelax = nrelax);

  foreach_face()
    ubf.x[] = -alpha.x[]*face_gradient_x (psi, 0);

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
psi[top]    = neumann (0);
psi[bottom] = neumann (0);
#  endif
#  if dimension > 2
psi[front]  = neumann (0);
psi[back]   = neumann (0);
#  endif
#endif // !AXI
