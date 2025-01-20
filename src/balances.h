#define BALANCES

extern int maxlevel;
extern scalar f, porosity;

struct MassBalances 
{
  //total mass variables
  double tot_sol_mass_start;
  double tot_gas_mass_start;
  double tot_sol_mass;
  double tot_gas_mass;
  double tot_mass_start;
  double tot_mass;
  double tot_gas_mass_bdnow;
  double tot_gas_mass_bd;

  //gas mass variables
  double* gas_massbdnow;
  double* gas_massbd;
  double* gas_mass_start;

  //file variables
  FILE* fb;
  char name[80];
  bool boundaries;
  int print_iter;
};

struct MassBalances mb = {0};

static double face_gradient_bid (Point point, scalar Y, int bid) {
  double grad = 0.;
  switch (bid) {
    case 0: grad = (Y[1,0]  - Y[])/Delta;  break;  // right
    case 1: grad = (Y[-1,0] - Y[])/Delta;  break;  // left
    case 2: grad = (Y[0,1]  - Y[])/Delta;  break;  // top
    case 3: grad = (Y[0,-1] - Y[])/Delta;  break;  // bottom
  }
  return grad;
}

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

static void diffusion_boundary (Point point, int bid) {
  foreach_elem (YGList_G, jj) {
    scalar YG = YGList_G[jj];
    double gradYG = face_gradient_bid (point, YG, bid);
    scalar Dmix2  = Dmix2List_G[jj];
    double Dmix2v = Dmix2[];
#if AXI
    mb.gas_massbdnow[jj] -= rhoG*Dmix2v*gradYG*Delta*dt*y;
#else
    mb.gas_massbdnow[jj] -= rhoG*Dmix2v*gradYG*Delta*dt;
#endif
  }
  foreach_elem (YGList_S, jj) {
    scalar YG = YGList_S[jj];
    double gradYG = face_gradient_bid (point, YG, bid);
    scalar Dmix2  = Dmix2List_S[jj];
    double Dmix2v = Dmix2[];
#if AXI
    mb.gas_massbdnow[jj] -= rhoG*Dmix2v*gradYG*Delta*dt*y;
#else
    mb.gas_massbdnow[jj] -= rhoG*Dmix2v*gradYG*Delta*dt;
#endif
  }
}

scalar U[];

static void advection_boundary (Point point, int bid) {
  foreach_elem (YGList_G, jj) {
    scalar YG = YGList_G[jj];
    double fluxYG = face_flux_bid (point, YG, bid);
    mb.gas_massbdnow[jj] += rhoG*fluxYG*dt;
  }
  foreach_elem (YGList_S, jj) {
    scalar YG = YGList_S[jj];
    double fluxYG = face_flux_bid (point, YG, bid);
    mb.gas_massbdnow[jj] += rhoG*fluxYG*dt;
  }

  /**
  Linear interpolation approximation for the calculation of the face
  volume fraction. */

  //double ff = face_value_bid (point, f, bid);
  //mb.totmass2bdnow += rhoG*face_flux_bid (point, U, bid, unity=true)*(1. - ff)*dt;
  mb.tot_gas_mass_bdnow += rhoG*face_flux_bid (point, U, bid, unity=true)*dt;
}

static void write_balances(void) {
  //print total mass
  fprintf(mb.fb, "%-12.6f %-12.6f %-12.6f %-12.6f ", t,
                                        mb.tot_sol_mass / mb.tot_sol_mass_start,
                                        (mb.tot_gas_mass - mb.tot_gas_mass_start) / mb.tot_sol_mass_start,
                                        mb.tot_mass / mb.tot_mass_start);

  //print individual gas species mass
  foreach_elem (YGList_G, jj)
    fprintf(mb.fb, "%-12.6f", (mb.gas_mass[jj] - mb.gas_mass_start[jj]) / mb.tot_sol_mass_start);
  
  //print individual solid species mass
  foreach_elem (YSList, jj)
    fprintf(mb.fb, "%-12.6f", (mb.sol_mass[jj]) / mb.tot_sol_mass_start);
  fprintf(mb.fb, "\n");
  fflush(mb.fb);
}

static void compute_balances(void) {
  mb.tot_sol_mass = 0.;
  mb.tot_gas_mass = 0.;
  mb.tot_mass = 0.;
  mb.tot_gas_mass_bdnow = 0.;

  foreach (serial) {
    mb.tot_sol_mass += (f[] - porosity[])*rhoS*dv(); //(f-ef)
    mb.tot_gas_mass += (1-f[] + porosity[])*rhoG*dv(); //(1-f + ef) 
  }

  if (mb.boundaries) {
    for (int b=0; b<nboundary; b++) {
      foreach_boundary (b) {
        diffusion_boundary (point, b);
        advection_boundary (point, b);
      }
    }
  }
  mb.tot_gas_mass_bd += mb.tot_gas_mass_bdnow;
  mb.tot_gas_mass += mb.tot_gas_mass_bd;
  mb.tot_mass = mb.tot_sol_mass + mb.tot_gas_mass;
}

event defualts (i = 0) {
  mb.tot_gas_mass_bdnow = 0.;
  mb.tot_gas_mass_bd = 0.;
  mb.boundaries = true;

  mb.gas_mass_start = (double*) malloc(NGS*sizeof(double));
  mb.gas_massbd =     (double*) malloc(NGS*sizeof(double));
  mb.gas_massbdnow =  (double*) malloc(NGS*sizeof(double));

  foreach()
    U[] = 1.;
}

event init(i = 0) {
  sprintf(mb.name, "balances-%d", maxlevel);
  mb.fb = fopen(mb.name, "w");
  //total mass header
  fprintf(mb.fb, "%-12s %-12s %-12s %-12s ", "t(1)", "tot_sol_mass(2)", "tot_gas_mass(3)", "tot_mass(4)");
  int counter = 5;
  //individual gas species mass header
  foreach_elem (YGList_G, jj) {
    char buffer[20];
    sprintf(buffer, "mG_%d(%d)", jj, counter++);
    fprintf(mb.fb, "%-12s ","mG", jj);
  }

  //solid species mass header
  foreach_elem (YSList, jj) {
    char buffer[20];
    sprintf(buffer, "mS_%d(%d)", jj, counter++);
    fprintf(mb.fb, "%-12s ","mS", jj);
  }
  
  fprintf(mb.fb, "\n");
  fflush(mb.fb);
  mb.tot_sol_mass_start = 0.;
  mb.tot_gas_mass_start = 0.;
  mb.tot_mass_start = 0.;
  mb.print_iter = 10;

  foreach (serial) {
    mb.tot_sol_mass_start += (f[] - porosity[])*rhoS*dv();
    mb.tot_gas_mass_start += (1-f[] + porosity[])*rhoG*dv();
  }
  mb.tot_mass_start = mb.tot_sol_mass_start + mb.tot_gas_mass_start;
  compute_balances();
}

event cleanup(t = end) {
  fclose(mb.fb);
}

event balances(i++, last) {
  compute_balances();
  if (i % mb.print_iter == 0)
    write_balances();
}