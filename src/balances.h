#define BALANCES

extern int maxlevel;
extern scalar f, porosity;

struct MassBalances 
{
  double sol_mass_start;
  double gas_mass_start;
  double sol_mass;
  double gas_mass;
  double total_mass_start;
  double total_mass;
  double gas_mass_bdnow;
  double gas_mass_bd;
  FILE* fb;
  char name[80];
  double* gas_massbdnow;
  bool boundaries;
};

struct MassBalances mb = {0};

// static double face_gradient_bid (Point point, scalar Y, int bid) {
//   double grad = 0.;
//   switch (bid) {
//     case 0: grad = (Y[1,0]  - Y[])/Delta;  break;  // right
//     case 1: grad = (Y[-1,0] - Y[])/Delta;  break;  // left
//     case 2: grad = (Y[0,1]  - Y[])/Delta;  break;  // top
//     case 3: grad = (Y[0,-1] - Y[])/Delta;  break;  // bottom
//   }
//   return grad;
// }

static double face_value_bid (Point point, scalar Y, int bid) {
  double val = 0.;
  switch (bid) {
    case 0: val = 0.5*(Y[1,0]  + Y[]); break;  // right
    case 1: val = 0.5*(Y[-1,0] + Y[]); break;  // left
    case 2: val = 0.5*(Y[0,1]  + Y[]); break;  // top
    case 3: val = 0.5*(Y[0,-1] + Y[]); break;  // bottom
  }
  return val;
}

static double face_flux_bid (Point point, scalar Y, int bid, bool unity=false) {
  double flux = 0., faceval = (unity) ? 1. : face_value_bid (point, Y, bid);
  switch (bid) {
    case 0: flux = +faceval*uf.x[]*Delta; break;  // right
    case 1: flux = -faceval*uf.x[]*Delta; break;  // left
    case 2: flux = +faceval*uf.y[]*Delta; break;  // top
    case 3: flux = -faceval*uf.y[]*Delta; break;  // bottom
  }
  return flux;
}

// static void diffusion_boundary (Point point, int bid) {
//   foreach_elem (YGList_G, jj) {
//     scalar YG = YGList_G[jj];
//     double gradYG = face_gradient_bid (point, YG, bid);
// #ifdef VARPROP
//     scalar Dmix2  = Dmix2List_G[jj];
//     double Dmix2v = Dmix2[];
// #else
//     double Dmix2v = mb.inDmix2[jj];
// #endif
// #if AXI
//     mb.gas_massbdnow[jj] -= rho2v[]*Dmix2v*gradYG*Delta*dt*y;
// #else
//     mb.gas_massbdnow[jj] -= rho2v[]*Dmix2v*gradYG*Delta*dt;
// #endif
//   }
//   foreach_elem (YGList_S, jj) {
//     scalar YG = YGList_S[jj];
//     double gradYG = face_gradient_bid (point, YG, bid);
// #ifdef VARPROP
//     scalar Dmix2  = Dmix2List_S[jj];
//     double Dmix2v = Dmix2[];
// #else
//     double Dmix2v = mb.inDmix2[jj];
// #endif
// #if AXI
//     mb.gas_massbdnow[jj] -= rho2v[]*Dmix2v*gradYG*Delta*dt*y;
// #else
//     mb.gas_massbdnow[jj] -= rho2v[]*Dmix2v*gradYG*Delta*dt;
// #endif
//   }
// }

scalar U[];

static void advection_boundary (Point point, int bid) {
  // foreach_elem (YGList_G, jj) {
  //   scalar YG = YGList_G[jj];
  //   double fluxYG = face_flux_bid (point, YG, bid);
  //   mb.gas_massbdnow[jj] += rhoG*fluxYG*dt;
  // }
  // foreach_elem (YGList_S, jj) {
  //   scalar YG = YGList_S[jj];
  //   double fluxYG = face_flux_bid (point, YS, bid);
  //   mb.gas_massbdnow[jj] += rhoG*fluxYG*dt;
  // }

  /**
  Linear interpolation approximation for the calculation of the face
  volume fraction. */

  //double ff = face_value_bid (point, f, bid);
  //mb.totmass2bdnow += rhoG*face_flux_bid (point, U, bid, unity=true)*(1. - ff)*dt;
  mb.gas_mass_bdnow += rhoG*face_flux_bid (point, U, bid, unity=true)*dt;
}

static void write_balances(void) {
  fprintf(mb.fb, "%-12.6f %-12.6f %-12.6f %-12.6f\n", t,
                                        mb.sol_mass / mb.sol_mass_start,
                                        (mb.gas_mass - mb.gas_mass_start) / mb.sol_mass_start,
                                        mb.total_mass / mb.total_mass_start);
  fflush(mb.fb);
}

static void compute_balances(void) {
  mb.sol_mass = 0.;
  mb.gas_mass = 0.;
  mb.total_mass = 0.;
  mb.gas_mass_bdnow = 0.;

  foreach (serial) {
    mb.sol_mass += (f[] - porosity[])*rhoS*dv(); //(f-ef)
    mb.gas_mass += (1-f[] + porosity[])*rhoG*dv(); //(1-f + ef) 
  }

  if (mb.boundaries) {
    for (int b=0; b<nboundary; b++) {
      foreach_boundary (b) {
        //diffusion_boundary (point, b);
        advection_boundary (point, b);
      }
    }
  }
  mb.gas_mass_bd += mb.gas_mass_bdnow;
  mb.gas_mass += mb.gas_mass_bd;
  mb.total_mass = mb.sol_mass + mb.gas_mass;
}

event defualts (i = 0) {
  mb.gas_mass_bdnow = 0.;
  mb.gas_mass_bd = 0.;
  mb.boundaries = true;

  foreach()
    U[] = 1.;
}

event init(i = 0) {
  mb.fb = fopen(mb.name, "w");
  fprintf(mb.fb, "%-12s %-12s %-12s %-12s\n", "t", "sol_mass", "gas_mass", "total_mass");
  fflush(mb.fb);
  mb.sol_mass_start = 0.;
  mb.gas_mass_start = 0.;
  mb.total_mass_start = 0.;
  foreach (serial) {
    mb.sol_mass_start += (f[] - porosity[])*rhoS*dv();
    mb.gas_mass_start += (1-f[] + porosity[])*rhoG*dv();
  }
  mb.total_mass_start = mb.sol_mass_start + mb.gas_mass_start;
  compute_balances();
}

event cleanup(t = end) {
  fclose(mb.fb);
}

event balances(i++, last) {
  compute_balances();
  if (i % 10 == 0)
    write_balances();
}