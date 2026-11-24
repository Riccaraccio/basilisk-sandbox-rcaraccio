/**
# Export the gas cells of a snapshot in the format of `CHEMISTRY_DUMP`

    export-cells DUMP KINFOLDER OUT

reads the Basilisk snapshot `DUMP` and writes the gas cells that the gas
reactor integrates (the test of `gas_reactor_cell()` in `chemistry.h`) to
`OUT.cells` and `OUT.meta`. The files have the format of `CHEMISTRY_DUMP`
(`src/chemistry-cell.h`), so `test/drg-dump.c` and cuSMOKE's
`make_reference --cells` can read them.

`KINFOLDER` is the folder of the kinetic scheme of the run, for example
`$OPENSMOKE_INTERFACE/kinetics/skeletal/isobuthanol-solid/kinetics`. The
exporter reads the names of the gas species from it. The fields of the
snapshot are `<species>_G` and `TG`, in the VOF-tracer form: the exporter
divides them by `1 - f`.

The `.meta` file holds P = 101325 + p for each cell, and zero for the flag
of the gate, the CPU time and the number of evaluations. A snapshot does not
hold these values.

Build:

    qcc -O2 -disable-dimensions -I$OPENSMOKE_INTERFACE/src export-cells.c \
      -L$OPENSMOKE_INTERFACE/build -lopensmoke -lm -o export-cells
*/

#include "utils.h"
#include "OpenSMOKE_Interface.h"
#include <stdint.h>

#define F_ERR 1e-10
#define PREF 101325.

int main (int argc, char * argv[]) {
  if (argc != 4) {
    fprintf (stderr, "usage: export-cells DUMP KINFOLDER OUT\n");
    return 1;
  }

  OpenSMOKE_Init ();
  OpenSMOKE_ReadKinetics (argv[2]);
  int ngs = OpenSMOKE_NumberOfSpecies ();

  init_grid (1);
  if (!restore (file = argv[1], list = all)) {
    fprintf (stderr, "export-cells: cannot restore %s\n", argv[1]);
    return 1;
  }

  scalar * YG = NULL;
  for (int jj = 0; jj < ngs; jj++) {
    char name[80];
    snprintf (name, sizeof (name), "%s_G", OpenSMOKE_NamesOfSpecies (jj));
    scalar s = lookup_field (name);
    if (s.i < 0) {
      fprintf (stderr, "export-cells: no field %s in the snapshot\n", name);
      return 1;
    }
    YG = list_append (YG, s);
  }
  scalar TG = lookup_field ("TG"), f = lookup_field ("f"),
    pf = lookup_field ("p");
  if (TG.i < 0) {
    fprintf (stderr, "export-cells: no field TG in the snapshot\n");
    return 1;
  }

  char name[1024];
  snprintf (name, sizeof (name), "%s.cells", argv[3]);
  FILE * fc = fopen (name, "wb");
  snprintf (name, sizeof (name), "%s.meta", argv[3]);
  FILE * fm = fopen (name, "wb");
  if (!fc || !fm) {
    fprintf (stderr, "export-cells: cannot open the output files\n");
    return 1;
  }

  /**
  The header: the number of cells (written again at the end), the number
  of species and the names in 32 bytes each. */

  int64_t ncells = 0;
  int32_t nsp = ngs;
  fwrite (&ncells, sizeof (int64_t), 1, fc);
  fwrite (&nsp, sizeof (int32_t), 1, fc);
  for (int jj = 0; jj < ngs; jj++) {
    char sp[32] = {0};
    strncpy (sp, OpenSMOKE_NamesOfSpecies (jj), 31);
    fwrite (sp, 1, 32, fc);
  }

  double row[4 + ngs];
  foreach (serial) {
    double ff = f.i >= 0 ? f[] : 0.;
    if (!(ff < 1. - F_ERR))
      continue;
    double temperature = TG[]/(1. - ff);
    if (!(temperature > 273.) || !(temperature < 3500.))
      continue;
    double ysum = 0.;
    int jj = 0;
    for (scalar s in YG) {
      row[4 + jj] = s[]/(1. - ff);
      ysum += s[];
      jj++;
    }
    if (!(ysum > 0.))
      continue;

    row[0] = x, row[1] = y;
    row[2] = (1. - ff)*sq(Delta);
    row[3] = temperature;
    fwrite (row, sizeof (double), 4 + ngs, fc);

    double meta[4] = {PREF + (pf.i >= 0 ? pf[] : 0.), 0., 0., 0.};
    fwrite (meta, sizeof (double), 4, fm);
    ncells++;
  }

  fseek (fc, 0, SEEK_SET);
  fwrite (&ncells, sizeof (int64_t), 1, fc);
  fclose (fc), fclose (fm);
  free (YG);

  fprintf (stderr, "%s: t = %g, %ld gas cells, %d species\n",
           argv[1], t, (long) ncells, ngs);
  return 0;
}
