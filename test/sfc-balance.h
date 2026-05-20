
/**
## New partition policy

Given a leaf's Z-order index, the global leaf count, the
number of processes and the cumulative-weight array, return the
process id
*/

static int new_balanced_pid (long index, long nt, int nproc, double* cf) {
  double il = cf[nt - 1]/nproc; // ideal load per process
  int idx = (int) floor ((cf[index] - cf[0])/il);
  return min (idx, nproc - 1);
}

/**
## Cumulative weight

`cf[i]` is the cumulative sum of `w` over leaves, in the same
order they are visited by the partitioning loop. */

void compute_cf (double* cf, long nt, scalar w) {
  long idx = 0;
  double acc = 0.;
  foreach_leaf() {
    acc += w[];
    cf[idx++] = acc;
  }
}

/**
## Per-process cell count and load 
Performance monitor function that returns the number of cells 
and the load in each processor*/

void balance_score (long* counter, double* load, scalar w) {
  for (int ne = 0; ne < npe(); ne++) {
    counter[ne] = 0;
    load[ne] = 0.;
  }

  foreach() {
    counter[pid()] += 1;
    load[pid()] += w[];
  }

  if (pid() == 0) {
    MPI_Reduce(MPI_IN_PLACE, counter, npe(), MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, load, npe(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  } else {
    MPI_Reduce(counter, NULL, npe(), MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(load, NULL, npe(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }
}

/**
## Static partitioning using the new policy

Copy of `mpi_partitioning()` with `balanced_pid()` replaced by
`new_balanced_pid()`. */

trace
void new_mpi_partitioning (scalar w)
{
  prof_start ("new_mpi_partitioning");

  long nt = 0;
  foreach (serial) nt++;

  //
  double tot_load = 0.;
  //foreach (reduction (+:tot_load)) //fixme why does this not work
  foreach (serial)
    tot_load += w[];
  //fprintf (stdout, "tot_load: %g\n", tot_load);

  double target_load = tot_load/npe();
  //fprintf (stdout, "#target_load: %g\n", target_load);
  int current_rank = 0;
  double running_load = 0.;
  //

  //long i = 0;
  tree->dirty = true;
  foreach_cell_post (is_active (cell))
    if (is_active (cell)) {
      if (is_leaf (cell)) {

        cell.pid = current_rank;
        running_load += w[];
        //fprintf (stdout, "%ld %g\n", i++, running_load);

        if (running_load >= target_load && current_rank < npe() - 1) {
          current_rank += 1;
          running_load = 0.;
        }

        //cell.pid = new_balanced_pid (i++, nt, npe(), cf);
        if (cell.neighbors > 0) {
          int pid = cell.pid;
          foreach_child() cell.pid = pid;
        }
        if (!is_local (cell)) cell.flags &= ~active;
      }
      else {
        cell.pid = child(0).pid;
        bool inactive = true;
        foreach_child() if (is_active (cell)) { inactive = false; break; }
        if (inactive) cell.flags &= ~active;
      }
    }
  flag_border_cells();
  prof_stop();
  mpi_boundary_update_buffers();
}

double parallel_efficency (double* load) {
  //compute parallel efficency given the load array
  double tot_load = 0., max_load = -1;
  for (int ne = 0; ne < npe(); ne++) {
    tot_load += load[ne];
    if (load[ne] > max_load)
      max_load = load[ne];
  }
  return (tot_load/npe())/max_load;
}
