/*-----------------------------------------------
  ____  _       ____            _ _ _     _    
 |  _ \(_)     |  _ \          (_) (_)   | |   
 | |_) |_  ___ | |_) | __ _ ___ _| |_ ___| | __
 |  _ <| |/ _ \|  _ < / _` / __| | | / __| |/ /
 | |_) | | (_) | |_) | (_| \__ \ | | \__ \   < 
 |____/|_|\___/|____/ \__,_|___/_|_|_|___/_|\_\

-----------------------------------------------*/

#include "myutils.h"
#define DIFFUSIVE 1
#define ufext uf
#define VARPROP
//#define NO_1D_COMPRESSION 1
//#define SOLVE_TEMPERATURE
//#define TURN_OFF_REACTIONS

#include "navier-stokes/centered-evaporation.h"
#include "two-phase.h"
#include "evaporation.h"
#include "memoryallocation.h"
#include "reactions.h"
//#include "osdiffusion.h"
#include "view.h"

#ifndef
# define circle(x, y, R) (sq(R*L0) - sq(x - 0.5*L0) - sq(y - 0.5*L0))
#endif

int MAXLEVEL = 6;
double R0 = 0.2;
double tEnd = 100 ;

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);

u.n[bottom] = neumann (0.);
u.t[bottom] = neumann (0.);
p[bottom] = dirichlet (0.);

u.n[left] = neumann (0.);
u.t[left] = neumann (0.);
p[left] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);

int main(){
  L0 = 1;
  DT = 0.1;
  rho1 = 845., rho2 = 1.25;
  mu1 = 0.4, mu2 = 1.e-5;
  cp1 = 1.76 , cp2 = 1;
  init_grid(1 << MAXLEVEL);
  kinfolder = "biomass-kinetics";
  f.tracers = {epsilon};
  run();
}

event init (i = 0){
  fraction (f, circle(x, y, R0)); 

  // Definition of non 0 species
  //GAS  
  //gas_start[OpenSMOKE_IndexOfSpecies ("O2")] = 0.21;
  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1;
  //SOLID 
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("CELL")-NGS] = 1;

  // Definiton of the temperature field
  foreach(){
    T[] = 650;
    epsilon[] = 0.6*f[];
  }
  Pref = 101325;
}

event movie (t += 1; t <= tEnd) {

  clear();
  view (tx = -0.5, ty = -0.5);
  draw_vof ("f");
  squares ("epsilon", min = 0, max = 0.6);  
  //labels("CELL");
  save ("movie.mp4");

}
