/**
# Negative values from the corrective flux of the coupled solve

The coupled solve linearizes the corrective term as
$$
Y_i U \approx [Y_i^{n+1}]_f U^n + [Y_i^n]_f (U^{n+1} - U^n)
$$
The second term takes the upwind value for the sign of $U^n$. Where $U$
changes in the step, it can take mass out of a cell with $Y_i = 0$.

The test uses a uniform grid, thus no change of level. A front of the
molecular weight (a light species `Ya` against a heavy carrier `Yb`) gives a
large corrective flux, as at a flame. A trace species `Yc` has a step or a
thin peak at a distance `s` from the front. The test gives the smallest
value of each species after one step, for a range of cell Fourier numbers
$Fo = D_a \Delta t/\Delta^2$, in the Fick form and in the molar form. */

#include "grid/multigrid.h"
#include "run.h"
#include "coupled-diffusion.h"

#define N 128
scalar Ya[], Yb[], Yc[], Da[], Db[], Dc[], thetac[];
face vector alpha[];
double MWs[3] = {2., 28., 162.};

/**
`w` is the width of the front in cells, `s` the position of the trace
species from the front in cells, `peak` a peak of one cell instead of a
step. */

static void initial (double w, double s, bool peak)
{
  double d = 1./N;
  foreach() {
    double xi = (x - 0.5)/d;
    Ya[] = 0.45 + 0.4*tanh (xi/w);
    if (peak)
      Yc[] = (xi > s - 0.5 && xi < s + 0.5) ? 1e-3 : 0.;
    else
      Yc[] = xi < s ? 1e-3 : 0.;
    Yb[] = 1. - Ya[] - Yc[];
    Da[] = 1., Db[] = 0.2, Dc[] = 0.1;
    thetac[] = cm[];
  }
  foreach_face()
    alpha.x[] = fm.x[];
}

int main()
{
  init_grid (N);
  scalar * YL = {Ya, Yb, Yc}, * DL = {Da, Db, Dc};
  double Fos[] = {0.25, 1., 4., 16., 64.};
  double ws[] = {1., 4.};
  double ss[] = {-4., -1., 0., 1., 4.};

  for (int molar = 0; molar <= 1; molar++)
    for (int peak = 0; peak <= 1; peak++)
      for (int iw = 0; iw < 2; iw++)
        for (int is = 0; is < 5; is++) {
          fprintf (stderr, "%s %s w=%g s=%+g:", molar ? "molar" : "Fick ",
                   peak ? "peak" : "step", ws[iw], ss[is]);
          for (int iF = 0; iF < 5; iF++) {
            double dt = Fos[iF]*sq(1./N);
            initial (ws[iw], ss[is], peak);
            mgstats st = coupled_diffusion (YL, dt, DL, alpha, thetac,
                                            MWs = molar ? MWs : NULL,
                                            tolerance = 1e-13/dt);
            double m[3] = {HUGE, HUGE, HUGE}, es = 0.;
            foreach (reduction(min:m[:3]) reduction(max:es)) {
              double sum = 0.;
              int k = 0;
              for (scalar Y in YL)
                m[k] = min (m[k], Y[]), sum += Y[], k++;
              es = max (es, fabs (sum - 1.));
            }
            fprintf (stderr, "  Fo %-4g c %+.1e a %+.0e b %+.0e%s",
                     Fos[iF], m[2]/1e-3, m[0], m[1],
                     st.resa > st.resb ? " DIV" :
                     st.resa > 1e-13/dt ? " NC" : "");
          }
          fprintf (stderr, "\n");
        }
}
