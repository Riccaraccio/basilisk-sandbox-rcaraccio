#define NO_ADVECTION_DIV 1
#define SOLVE_TEMPERATURE 1
#define MULTICOMPONENT 1

#include "axi.h"
#include "navier-stokes/centered-phasechange.h"
#include "constant-properties.h"
#include "two-phase.h"
#include "shrinking.h"
#include "multicomponent-varprop.h"
#include "darcy.h"
//#include "flame.h"

const double Uin = 0.13; //inlet velocity
u.n[left]    = dirichlet (Uin);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
psi[left]    = dirichlet (0.);

psi[top]     = dirichlet (0.);

u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
psi[right]    = neumann (0.);

double tend = 150;
int maxlevel = 9, minlevel = 2;
double solid_mass0 = 0.;
double D0 = 8e-3;

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

int main() {

  TS0 = 300.; TG0 = 1123.;
  rhoS = 1550;
  eps0 = 0.2;

  rhoG    = 0.31;      // air ~1100 K, 1 atm
  muG     = 4.5e-5;
  lambdaG = 0.08;      cpG = 1200.;
  lambdaS = 0.2;       cpS = 1500.;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_SWELLING;

  DT = 1e-2;

  kinfolder = "biomass/dummy-solid";
  shift_prod = true;

  L0 = 20*D0;
  origin (-L0/2, 0);

  emissivity = emissivity_constant;

  init_grid(1 << min (maxlevel, 8));
  refine (circle (x, y, 4.*D0) > 0. && level < maxlevel);

  run();
}

double r0;
event init (i= 0) {
  scalar f0[];

  fraction (f0, circle(x, y, 0.5 * D0));

  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;
  //gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 0.765;
  //gas_start[OpenSMOKE_IndexOfSpecies ("O2")] = 0.235;

  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")] = 1.;

  foreach()
    porosity[] = eps0*f0[];

  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0))
    solid_mass0 += f0[]*(1. - eps0)*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f

  TG[left] = dirichlet (TG0);
  TG[top] = dirichlet (TG0);

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) {
      //YG[left] = dirichlet (0.765);
      //YG[top] = dirichlet (0.765);
      YG[left] = dirichlet (1.);
      YG[top] = dirichlet (1.);
    } 
    //else if (jj == OpenSMOKE_IndexOfSpecies ("O2")) {
    //  YG[left] = dirichlet (0.235);
    //  YG[top] = dirichlet (0.235);
    //}
    else {
      YG[left] = dirichlet (0.);
      YG[top] = dirichlet (0.);
    }
  }

  if (restore (file = "last-snapshot", list = all)) {
    fprintf (stderr, "Restart file found!\n");
    restarted = true;
  } else {
    fprintf (stderr, "No restart file found, starting from scratch!\n");

    foreach() {
      f[] = f0[];
      porosity[] = eps0*f[];
    }
  }
}

// Calculates the H2O-density path-averaged temperature
double T_H2O_weigthed_average (double x_interp, int n_samples = 200, const double length = L0/2) {
  scalar YH2O = YGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];
  //scalar XH2O = XGList_G[OpenSMOKE_IndexOfSpecies ("H2O")];

  double numerator = 0., denominator = 0.;
  coord pos, box[2] = {{x_interp, 0.}, {x_interp, length}}, nn = {1, n_samples};
  foreach_region (pos, box, nn, reduction(+:numerator) reduction(+:denominator)) {
    double yH2O_local = interpolate_linear (point, YH2O, pos.x, pos.y, pos.z);
    numerator   += yH2O_local;
    denominator += yH2O_local / interpolate_linear (point, T, pos.x, pos.y, pos.z);
  }

  if (denominator <= 0) // avoid division by 0
    return TG0;

  return numerator/denominator;
}

event output (t += 0.01) {

  char name[80];
  sprintf(name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, restarted ? "a" : "w");
  if (fp == NULL) {
    fprintf (stderr, "Error opening OutputData\n");
    exit(1);
  }

  // Interpolate Water vapor mass fraction
  double Tavg[3], sample_points[3] = {D0/2 + 2e-3, D0/2 + 4e-3, D0/2 + 11e-3};

  const double L_flame_exp = 20e-3/2;

  for (int ii = 0; ii < 3; ii++)
      Tavg[ii] = T_H2O_weigthed_average (sample_points[ii]);

  if (i == 0)
    fprintf (fp, "#t(1), Ms/Ms0(2), Tmax(3), Tavg_1mm(4), Tavg_2mm(5), Tavg_4mm(6)\n");

  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[] - porosity[])*rhoS*dv();

  fprintf (fp, "%g %g %g %g %g %g\n", t, solid_mass/solid_mass0, statsf(T).max,
                                      Tavg[0], Tavg[1], Tavg[2]);

  fflush(fp);
}

/**
Diagnostics for the phase-change expansion field, sampled every timestep.

The projection enforces div(uf) = -gas_source, so we log the two sides of that
identity separately:
  Qsrc = int gas_source dV   the source computed by the chemistry
  Qdiv = int div(uf) dV      the expansion the velocity field actually carries
and 'resmax', the largest pointwise violation. Note gas_source already carries
cm[], as does the discrete div(uf) below, so both integrate with sq(Delta).

'ur' probes the radial velocity on a 45 degree ray just outside the particle:
this is the quantity the velocity vectors show.
*/

event probe_expansion (i++) {
  double Qsrc = 0., Qdiv = 0., resmax = 0.;

  foreach (reduction(+:Qsrc) reduction(+:Qdiv) reduction(max:resmax)) {
    double d = 0.;
    foreach_dimension()
      d += uf.x[1] - uf.x[];
    d /= Delta;                        // = cm*div(u), same weighting as gas_source

    Qsrc += gas_source[]*sq(Delta);
    Qdiv += d*sq(Delta);
    resmax = max (resmax, fabs (d + gas_source[]));
  }

  stats so = statsf (omega);

  double ur[3];
  for (int k = 0; k < 3; k++) {
    double r = 0.5*D0 + (k + 1)*0.5e-3, c = cos(pi/4.), s = sin(pi/4.);
    ur[k] = interpolate (u.x, r*c, r*s)*c + interpolate (u.y, r*c, r*s)*s;
  }

  /**
  Particle-side temperature. statsf(T).max is useless here: it sits on TG0
  because the particle starts at TS0 < TG0 and is always the coldest region,
  so a maximum only ever reports the inlet boundary condition.

  We use T[] (assigned once per step as TS[]+TG[], both in tracer form, so it
  is the volume-weighted mixture temperature) together with f[]. Both are
  unambiguous wherever this event runs, unlike TS[]/TG[] whose normalisation
  depends on the position within the timestep.

    Tcore  coldest point, i.e. the centre of the particle
    Tbulk  f-weighted mean over the solid, the bulk particle temperature
    Tsurf  mean over interface cells: where the Arrhenius feedback bites first
  */

  double Tcore = statsf(T).min;
  double Tbulk = 0., fvol = 0., Tsurf = 0., nsurf = 0.;

  foreach (reduction(+:Tbulk) reduction(+:fvol)
           reduction(+:Tsurf) reduction(+:nsurf)) {
    if (f[] > F_ERR) {
      Tbulk += T[]*f[]*dv();
      fvol  += f[]*dv();
    }
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      Tsurf += T[];
      nsurf += 1.;
    }
  }
  Tbulk = fvol   > 0. ? Tbulk/fvol  : 0.;
  Tsurf = nsurf  > 0. ? Tsurf/nsurf : 0.;

  /**
  Everything above is collective (reductions / interpolate), so every rank
  holds the same values. Only rank 0 writes: basilisk wraps popen for MPI but
  not fopen, so an unguarded fprintf would emit one copy of each line per rank.
  */

  if (pid() == 0) {
    static FILE * fe = fopen ("expansion.dat", restarted ? "a" : "w");
    if (fe == NULL) {
      fprintf (stderr, "Error opening expansion.dat\n");
      exit(1);
    }
    if (i == 0)
      fprintf (fe, "#t(1) dt(2) Qsrc(3) Qdiv(4) resmax(5) omega_min(6) omega_max(7)"
                   " ur_0.5mm(8) ur_1mm(9) ur_1.5mm(10) mgp_i(11) mgp_resa(12)"
                   " mgpsf_i(13) Tcore(14) Tbulk(15) Tsurf(16)\n");

    fprintf (fe, "%g %g %g %g %g %g %g %g %g %g %d %g %d %g %g %g\n",
             t, dt, Qsrc, Qdiv, resmax, so.min, so.max,
             ur[0], ur[1], ur[2], mgp.i, mgp.resa, mgpsf.i,
             Tcore, Tbulk, Tsurf);
    fflush (fe);
  }
}

/**
## Angular profile along the particle surface

For a sphere blowing uniformly the *normal* velocity at the surface is the same
at every angle,

  un_pred = Q/(4 pi R^2),   Q = 2 pi Qdiv the total volumetric production,

which reduces to Qdiv/(2 R^2). |u| is *not* uniform even then, because the free
stream adds a tangential component that varies with angle. So 'un' answers "is
the blowing uniform?" while 'umag' is what the |u| screenshots show; logging
both separates the anomaly from the free-stream confound.

theta runs from the downstream pole (0 deg) through the equator (90) to the
upstream pole (180). Along each ray we locate the interface, sample the gas
just outside it, and take the *peak* of omega inside it -- peak rather than a
fixed depth, so the measurement follows the reaction front instead of sliding
off it as the front recedes. 'r_front' vs theta is then the front shape, and
'T_front' the temperature driving the local Arrhenius rate.

Every interpolate() here is collective and called an identical number of times
on every rank (the branch conditions depend only on reduced values), so the
arrays agree across ranks and only rank 0 writes.
*/

#define NANG 24        // angular samples over [0,pi]
#define NRAD 32        // radial samples along each ray

double profile_offset = 0.15;  // gas-side sampling offset, in units of R

event angular_profile (t += 0.01) {
  const double R = 0.5*D0, dr = profile_offset*R;

  double th[NANG], ri[NANG], un[NANG], um[NANG], om[NANG], rf[NANG], Tf[NANG];

  for (int k = 0; k < NANG; k++) {
    double theta = (k + 0.5)*pi/NANG, c = cos(theta), s = sin(theta);

    /**
    Locate the interface along the ray. Marching outward to the first f < 0.5
    keeps this correct if the interface ever starts moving again.
    */

    double rint = R;
    for (int j = 0; j < NRAD; j++) {
      double rr = (j + 0.5)*1.5*R/NRAD;
      double ff = interpolate (f, rr*c, rr*s);
      if (ff != nodata && ff < 0.5) { rint = rr; break; }
    }

    /**
    Peak reaction rate along the ray, and the temperature where it peaks.
    */

    double ommax = 0., rfront = 0.;
    for (int j = 0; j < NRAD; j++) {
      double rr = (j + 0.5)*rint/NRAD;
      double o = interpolate (omega, rr*c, rr*s);
      if (o != nodata && fabs(o) > fabs(ommax)) {
        ommax = o; rfront = rr;
      }
    }

    /**
    Gas side, just outside the interface.
    */

    double ux = interpolate (u.x, (rint + dr)*c, (rint + dr)*s);
    double uy = interpolate (u.y, (rint + dr)*c, (rint + dr)*s);

    th[k] = theta*180./pi;
    ri[k] = rint;
    un[k] = ux*c + uy*s;                 // outward normal component = blowing
    um[k] = sqrt (sq(ux) + sq(uy));
    om[k] = ommax;
    rf[k] = rfront;
    Tf[k] = interpolate (T, rfront*c, rfront*s);
  }

  /**
  Total production, for the uniform-blowing reference. Same weighting as in
  probe_expansion: the discrete divergence carries cm[], so it integrates
  with sq(Delta).
  */

  double Qdiv = 0.;
  foreach (reduction(+:Qdiv)) {
    double d = 0.;
    foreach_dimension()
      d += uf.x[1] - uf.x[];
    Qdiv += d*Delta;
  }
  double un_pred = Qdiv/(2.*sq(R));

  if (pid() == 0) {
    static FILE * fa = fopen ("angular.dat", restarted ? "a" : "w");
    if (fa == NULL) {
      fprintf (stderr, "Error opening angular.dat\n");
      exit(1);
    }
    if (i == 0)
      fprintf (fa, "#t(1) theta_deg(2) r_int(3) un(4) umag(5) omega_max(6)"
                   " r_front(7) T_front(8) un_pred(9)\n");

    for (int k = 0; k < NANG; k++)
      fprintf (fa, "%g %g %g %g %g %g %g %g %g\n",
               t, th[k], ri[k], un[k], um[k], om[k], rf[k], Tf[k], un_pred);
    fflush (fa);
  }
}

#if TREE
event adapt (i++) {
  //scalar oxidiser = YGList_G[OpenSMOKE_IndexOfSpecies ("O2")];
  scalar oxidiser = YGList_G[OpenSMOKE_IndexOfSpecies ("N2")];

  //scalar zdiff[];
  //foreach()
  //  zdiff[] = zmix[] - zsto[];

  //adapt_wavelet_leave_interface ({T, oxidiser, zdiff}, {f},
  //  (double[]){5e0, 1.e-2, 1.e-2}, maxlevel, minlevel, 2);
  
  adapt_wavelet_leave_interface ({T, oxidiser}, {f}, (double[]){5e0, 1.e-2}, maxlevel, minlevel, 2);

  // Unrefine for outflow condition
  unrefine (x > L0*0.4);
}
#endif

//event movie (t += 0.1) {
//  clear();
//  view (theta=0, phi=0, psi=-pi/2., width = 1080, height = 1080);
//  squares ("T", min = 300, max = 2000, spread = -1, linear = true);
//  isoline ("T", val = statsf(T).max);
//  isoline ("zmix - zsto", lw = 1.5, lc = {1., 1., 1.});
//  draw_vof ("f", lw = 1.5);
//  mirror ({0, 1}) {
//    squares ("O2_G + O2_S", min = 0., max = 0.235, spread = -1, linear = true);
//    isoline ("zmix - zsto", lw = 1.5, lc = {1., 0., 0.});
//    draw_vof ("f", lw = 1.5);
//  }
//  save ("movie.mp4");
//}

event dump (t = 1; t += 1) {
  dump("last-snapshot");
}

event stop (t = tend);

/**
~~~gnuplot
~~~
**/
