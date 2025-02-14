#ifndef DARCY
# define DARCY_INDEX 0
#endif

#define POROUS_MEDIA 1

#ifndef POROUS_MEDIA 
  #include "embed.h"
#endif
#include "navier-stokes/centered.h"
#ifdef POROUS_MEDIA
  #include "two-phase.h"
  #include "darcy.h"
#endif
#include "view.h"
#include "adapt_wavelet_leave_interface.h"

////INPUTS//////////////
int maxlevel = 10; 
double Re = 20;
double R0 = 0.5;
double U0 = 1.;
double epsi0 = 0.3;
double side_length = 40;
double tend = 150.;
////////////////////////

// BOUNDARIES: left and right open, top and bottom no-slip
u.n[left] = dirichlet (U0);
u.t[left] = dirichlet (0.);
p[left] = neumann (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);

#ifndef POROUS_MEDIA
  u.n[embed] = fabs(y) > R0 ? neumann(0.) : dirichlet(0.);
  u.t[embed] = fabs(y) > R0 ? neumann(0.) : dirichlet(0.);
#endif

#ifdef POROUS_MEDIA
  scalar porosity[];
  double rhoG, muG;
#else
  face vector muv[];
#endif

int main() {
  double DaList[] = {1.00E-05, 5.00E-05, 1.00E-04, 5.00E-04, 1.00E-03, 2.50E-03, 5.00E-03};
  #ifndef POROUS_MEDIA
    mu = muv;
  #else
    mu1 = R0*2*U0/Re;
    mu2 = R0*2*U0/Re;
    rhoG = 1;
    muG = mu1;
    Da = DaList[DARCY_INDEX];
    fprintf(stderr, "Da = %g\n", Da);
  #endif

  L0 = side_length*R0*2;
  DT = 1e-3;
  origin(-L0/2, 0);
  init_grid (1 << maxlevel);
  run();
}

#define circle(x, y, R) (sq(R) - sq(x) - sq(y))

#ifndef POROUS_MEDIA
  event properties (i++){
    foreach_face()
      muv.x[] = fm.x[]*R0*2*U0/Re;
  }
#endif

event init (i = 0) {
  mask (y > L0/2 ? top : none);
  #ifdef POROUS_MEDIA
    fraction (f, circle(x, y, R0));

    foreach()
      porosity[] = f[]*epsi0;
  
    foreach() {
      u.x[] = f[] > 0 ? 0 : U0; 
    }

    foreach_face(x){
      double ff = face_value(f, 0);
      uf.x[] = ff > 0 ? 0 : U0;
    }
    #else
      solid (cs, fs, -circle(x, y, R0));
      foreach()
        u.x[] = cs[] > 1.- 1.e-10 ? U0 : 0;
    #endif
}

#ifdef POROUS_MEDIA
  // Avoid tranporting the interface: the solid is fixed
  static scalar * interfaces_save = NULL;

  event vof (i++) {
    interfaces_save = interfaces; 
    interfaces = NULL;
  }
  event tracer_advection (i++) {  
    interfaces = interfaces_save;
  }
#endif

# if   AMR_ACTIVE
event adapt (i++) {
#  ifndef POROUS_MEDIA  
  adapt_wavelet_leave_interface ({u.x, u.y}, {cs}, (double[]){1.e-3, 1.e-3}, maxlevel, 3, padding=1);
#  else
  adapt_wavelet_leave_interface ({u.x, u.y}, {f}, (double[]){1.e-3, 1.e-3}, maxlevel, 3, padding=1);
#  endif
}
# endif

scalar psi[];
scalar omega[];
// Movie event
event movie (t += 1) {
  foreach()
    omega[] = 0;

  vorticity (u, omega);
  psi[top] = dirichlet(0.);
  psi[bottom] = dirichlet(U0*L0);  

  poisson(psi, omega);

  view (fov = 6, tx = -0.65, ty = -0.5,
        width = 1600, height = 800);
  box ();
  #ifdef POROUS_MEDIA
    draw_vof ("f", fc = {0.5,0.5,0.5});
  #else
    draw_vof ("cs", fc = {0.5,0.5,0.5});
  #endif
  isoline (phi = "psi", n = 80, min=0.4 , max=0.6);
  save("streamlines.mp4");

  clear();

  view (fov = 6, tx = -0.65, ty = -0.5,
        width = 1600, height = 800);
  box ();
  #ifdef POROUS_MEDIA
    draw_vof ("f", fc = {0.5,0.5,0.5});
  #else
    draw_vof ("cs", fc = {0.5,0.5,0.5});
  #endif
  squares (color = "u.x", spread = -1, linear = true);
  save("velocity.mp4");
} 

//Print the current time on log file
event logfile (t += 1) {
  fprintf(stderr,"%g \n", t);
}

// Output data
void write_profile (void){
  char name[80];
  #ifdef POROUS_MEDIA
    sprintf (name, "WakeLength-Da-%d", DARCY_INDEX);
  #else
    sprintf (name, "WakeLength-Embed");
  #endif
  FILE * fp = fopen (name, "w");

  double step = L0/(1<<maxlevel);
  for (double x = 0; x<L0; x+=step){
    fprintf(fp, "%g %g\n", x, interpolate(u.x, x, 0));
  }

  fprintf (fp, "\n\n");
  fflush(fp);
}

event end (t = tend)
{
  printf("Reached final time\n");
  write_profile();
}

#if DUMP
int count = 0;
event snapshots (t += 1) {
  // we keep overwriting the last two snapshots
  if (count == 0) {
    dump ("snapshot-0");
    count++;
  } else {
    dump ("snapshot-1");
    count = 0;
  }
}
#endif