#include "vofToLs.h"
#include "velocity-potential.h"
#include "common-evaporation.h"

double eps0 = 0.5;
double rhoS = 100.;
double rhoG = 1.;
double muG = 1e-5;

extern scalar omega;

scalar porosity[];
scalar zeta[];
scalar levelset[];
scalar o[];

extern double TG0, TS0;
extern scalar TS, T;
vector gTS[];
scalar modg[];
double o_max = 0.;
scalar ls[];


/////////////SMA

#define WINDOW_SIZE 1000

typedef struct {
  double window[WINDOW_SIZE];
  int count;
  int index;
  double sum;
} MovingAverage;

// Initialize the moving average structure
void init_moving_average(MovingAverage *ma) {
    ma->count = 0;
    ma->index = 0;
    ma->sum = 0.0;
}

// Add a new value and compute the moving average
double add_value(MovingAverage *ma, double value) {
    if (ma->count < WINDOW_SIZE) {
        ma->count++;
    } else {
        // Remove the oldest value from the sum
        ma->sum -= ma->window[ma->index];
    }

    // Add the new value to the window and sum
    ma->window[ma->index] = value;
    ma->sum += value;

    // Move the index in a circular fashion
    ma->index = (ma->index + 1) % WINDOW_SIZE;

    return ma->sum / ma->count;
}

MovingAverage ma;
////////////

void tracer_gradients(scalar tr, scalar f, vector g) {
  // Compute gradient of tracer field tr using centered differences
  // tr: tracer field, ASSUMED TO NOT BE IN TRACER FORM
  // f: volume fraction field
  // g: gradient field

  foreach() {
    if (f[] > F_ERR) {
      foreach_dimension() {
        if (f[-1] > F_ERR && f[1] < F_ERR) {
          g.x[] = (tr[]-tr[-1])/Delta;
        } else if (f[-1] < F_ERR && f[-1] > F_ERR) {
          g.x[] = (tr[1]-tr[])/Delta;
        } else {
          g.x[] = (tr[1]-tr[-1])/(2.*Delta);
        }
      }
    }
  }
}


int compare_double(const void *a, const void *b) {
  return (*(double*)a > *(double*)b) - (*(double*)a < *(double*)b);
}

double findMax (double *arr, int n) {
  qsort(arr, n, sizeof(double), compare_double);
  int index = (int)ceil(0.95 * (n-1));
  if (index == n)
    return 0.;
  return arr[index];
}

double avgTop (double *arr, int n) {
  qsort(arr, n, sizeof(double), compare_double);
  int index = (int)ceil(0.95 * (n-1));
  if (index == n)
    return 0.;

  double sum = 0.;
  for (int i = index; i < n; i++)
    sum += arr[i];
  
  return sum/(n-index);
}

enum zeta_types {
  ZETA_SHRINK,
  ZETA_SWELLING,
  ZETA_SMOOTH,
  ZETA_SHARP,
  ZETA_LEVELSET,
  ZETA_REACTION,
  ZETA_GRADIENT,
  ZETA_LAST,
  ZETA_AVGTOP
};

enum zeta_types zeta_policy;

void set_zeta (enum zeta_types zeta_policy) {
  #if AXI
    double radius = pow (3.*statsf(f).sum, 1./3.);
  #else
    double radius = sqrt (statsf(f).sum/pi)*2;
  #endif

  switch (zeta_policy) {
    case ZETA_SHRINK:
      foreach()
        zeta[] = 0.5;//temp
      break;

    case ZETA_SWELLING:
      foreach()
        zeta[] = 0.;
      break;

    case ZETA_SMOOTH:
      foreach()
          //zeta[] = 1 / (1 + exp(32*radius - 40*sqrt(sq(x)+sq(y)+sq(z))));
          zeta[] = 1 / (1 + exp(-(sqrt(sq(x)+sq(y)+sq(z))-radius/2)/pow(radius, 4./3.)));
      break;

    case ZETA_SHARP:
      foreach()
        zeta[] = (sqrt(sq(x) + sq(y)) > radius*0.8) ? 1. : 0.;
      break;

    case ZETA_LEVELSET: {
      // vof_to_ls (f, levelset, imax=5);
      // double lmin = statsf(levelset).min;
      // if (fabs(lmin) > F_ERR)
      //   foreach() {
      //     zeta[] = 1 - levelset[]/statsf(levelset).min;
      //     zeta[] = clamp(zeta[], 0., 1.);
      //   }
      
      // vof_to_ls (f, levelset, imax=60);
      break;
    }

    case ZETA_REACTION: {
      foreach()
        o[] = omega[]*f[];

      o_max = statsf(o).max;

      foreach() {
        zeta[] = o_max > F_ERR ? omega[]/o_max : 0.;
        zeta[] = clamp(zeta[], 0., 1.);
      }
      break;
    }
    case ZETA_GRADIENT: {
      #ifdef SOLVE_TEMPERATURE
      foreach()
        TS[] = f[] > F_ERR ? TS[] / f[] : 0.;

      //gradients({TS}, {gTS});
      gradients({T}, {gTS});
      // tracer_gradients(TS, f, gTS);

      foreach() {
        modg[] = 0.;
        foreach_dimension()
            modg[] += sq(gTS.x[]);
        modg[] = sqrt(modg[]);
      }

      foreach() {
        
        // double Lc = (D0*0.5)*H0/(2*(D0*0.5+H0));
        double Lc = 0.5*D0;
        double Rc = modg[]/fabs(TG0-TS0)*Lc;
       
        zeta[] = f[] > F_ERR ? Rc: 0.;
        zeta[] = clamp(zeta[], 0., 1.);

        // zeta[] = f[] > F_ERR ? tanh(modg[]/TS[]*delta) : 0.;
        // zeta[] = f[] > F_ERR ? tanh(modg[]/TS[]*delta) : 0.;
        // zeta[] = TS[] > 0. ? 1 - 1/((modg[]/TS[])+1) : 0.;

        TS[] = f[] > F_ERR ? TS[]*f[] : 0.;
      }
      #endif
      break;
    } 
    case ZETA_LAST: {
      double* omega_arr = NULL;
      int n = 0;

      foreach(serial)
        if (f[] > F_ERR) {
          omega_arr = (double*) realloc (omega_arr, (n+1)*sizeof(double));
          omega_arr[n] = omega[]*f[];
          n++;
        }

      double omega_max = findMax(omega_arr, n);
      free (omega_arr); 

      foreach() {
        if (omega_max > F_ERR) {
          zeta[] = f[] > F_ERR ? omega[]/omega_max : 0.;
          zeta[] = clamp(zeta[], 0., 1.);
        } else {
          zeta[] = 0.;
        }
      }
      
      break;
    }
    case ZETA_AVGTOP: {
      //compute moving window average of the maximun
      //store 
      foreach()
        o[] = omega[]*f[];

      o_max = add_value(&ma, statsf(o).max);

      foreach() {
        zeta[] = f[] > F_ERR ? omega[]/o_max : 0.;
        zeta[] = clamp(zeta[], 0., 1.);
      }

      break;
    }
    default:
      fprintf (stderr, "Unknown Shrinking model\n");
      return;
  }
}

event defaults (i=0) {
  f.tracers = list_append (f.tracers, porosity);
  init_moving_average(&ma);
  set_zeta (zeta_policy);
}

event reset_sources (i++);

event chemistry (i++);

event phasechange (i++) {

  foreach() {
    f[] = clamp(f[], 0.,1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
    porosity[] = clamp(porosity[], 0., 1.);
    porosity[] = (porosity[] > F_ERR) ? porosity[] : 0.;
  }

  set_zeta (zeta_policy);

  face vector alpha_save[];
  foreach_face() {
    alpha_save.x[] = alpha.x[];
    alpha.x[] = fm.x[];
  }

  mgpsf = project_sv (ubf, psi, alpha, dt, mgpsf.nrelax); //TODO: should recive fm instead of alpha

  foreach_face()
    alpha.x[] = alpha_save.x[];

  foreach() {
    if (f[] > F_ERR) {
#ifdef VARPROP //*(f-ef)
      gasSource[] = -omega[]*(f[]-porosity[])*(1/rhoGv_S[] - 1/rhoSv[])*cm[];
#else
      gasSource[] = -omega[]*(f[]-porosity[])*(1/rhoG - 1/rhoS)*cm[];
#endif
    }
  }
}

face vector ufsave[];
event vof(i++) {
  foreach_face()
    ufsave.x[] = uf.x[];

  foreach_face()
    uf.x[] = ubf.x[];
}

event tracer_diffusion (i++) {
  foreach_face()
    uf.x[] = ufsave.x[];
}

event stability (i++, last) {
  face vector us[];
  foreach_face()
    us.x[] = ubf.x[];
  dt = dtnext (timestep (us, dtmax));
}

event cleanup (t = end) {
  free(f.tracers), f.tracers = NULL;
}
