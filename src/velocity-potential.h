/**
# Velocity potential solver for the solid phase velocity field
The calculation of the solid phase velocity field 'ubf' is performed through the
solution of a Poisson equation for the velocity potential 'psi'.
The source term of the Poisson equation is given by the product of the
reaction rate 'omega', the volume fraction field 'f' and the
shrinkage factor 'zeta', divided by the solid phase density 'rhoS'.
*/

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

scalar prod[];
extern scalar fS;

/**
## Projection method for the solid phase velocity field
We solve the Poisson equation for the velocity potential 'psi' and
then compute the solid phase velocity field 'ubf' as the negative gradient
of 'psi'.
*/

/**
The source term has a nonzero integral. If all the boundary conditions of
'psi' are Neumann conditions, the Poisson problem has no solution.
*psi_dirichlet_faces()* gives the number of boundary faces with a Dirichlet
condition on 'psi'. It calls the homogeneous condition of each boundary with a
field equal to 1. A Dirichlet condition gives -1 in the ghost cell. A Neumann
condition, a symmetry condition, and the axis give +1.
*psi_check_dirichlet()* stops the run if this number is zero. */

int psi_dirichlet_faces (scalar psi)
{
  scalar one[];
  foreach()
    one[] = 1.;
  int ndirichlet = 0;
  for (int b = 0; b < 2*dimension; b++) {
    int nb = 0;
    foreach_boundary (b, reduction(+:nb))
      if (psi.boundary_homogeneous[b] (point, neighborp(ig,jg,kg), one, NULL) < 0.)
        nb++;
    ndirichlet += nb;
  }
  return ndirichlet;
}

static void psi_check_dirichlet (scalar psi)
{
  if (psi_dirichlet_faces (psi) == 0) {
    fprintf (stderr, "velocity-potential.h: psi has no Dirichlet boundary "
             "condition. Set psi[<open boundary>] = dirichlet (0.) in the case.\n");
    exit (1);
  }
}

trace
mgstats project_sv (face vector ubf, scalar psi,
    (const) face vector alpha = unityf,
    int nrelax = 4)
{
  static bool psi_bc_checked = false;
  if (!psi_bc_checked) {
    psi_check_dirichlet (psi);
    psi_bc_checked = true;
  }

  /**
  We compute the source term for the Poisson equation
  */

  foreach()
    prod[] = omega[]*f[]*zeta[]*cm[]/rhoS;
 
  /**
  We optionally shift and/or diffuse the source term to improve stability
  when the heat exchange at the interface is high and the reaction rate
  is very localized.
  */

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

  mgstats mgp = poisson (psi, prod, alpha,
      tolerance = TOLERANCE_SOLID, nrelax = nrelax);

  foreach_face()
    ubf.x[] = -alpha.x[]*face_gradient_x (psi, 0);

  return mgp;
}

/**
## Boundary conditions

The default is a zero normal solid velocity on each boundary. The case must
set a Dirichlet condition on at least one open boundary (see
*psi_check_dirichlet()*). Do not copy the conditions of the pressure
(*neumann_pressure()*): they come from the Navier--Stokes relation for 'p',
not from the relation $\partial\psi/\partial n = -u_{b,n}$ for 'psi'.
*/

psi[right] = neumann (0);
psi[left]  = neumann (0);

#if AXI
ubf.n[bottom] = 0.;
ubf.t[bottom] = dirichlet(0);
psi[top]    = neumann (0);
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

event defaults (i=0) {
  psi.nodump = true;
  #if TREE
  ubf.x.refine = refine_face_solenoidal;
  #endif
}
