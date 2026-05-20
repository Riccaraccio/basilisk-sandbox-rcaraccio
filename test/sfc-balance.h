
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

Copy of `mpi_partitioning()` with the new partitioning method*/

trace
void new_mpi_partitioning ((const) scalar w = {-1}) {
  prof_start ("new_mpi_partitioning");

  long nt = 0;
  foreach (serial) nt++;
  
  //
  // if weigths are not provided, we set unity field
  if (w.i < 0)
    w[] = 1.;

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


trace
bool new_balance((const) scalar w = {-1}) {
  if (npe() == 1)
    return false;

  assert (sizeof(NewPid) == sizeof(double));

  check_flags();

  //if weigths are not set, we set the unity field 
  if (w.i < 0)
    w[] = 1.;

  long nl = 0, nt = 0;
  double tot_load = 0.;
  foreach_cell() {
    if (is_local(cell)) {
      nt++;
      if (is_leaf(cell)) {
        nl++;
        tot_load += w[]; // fixme: is it correct to count only leaves?
      }
    }
    if (is_leaf(cell))
      continue;
  }

  grid->n = grid->tn = nl;
  grid->maxdepth = depth();
  long nmin = nl, nmax = nl;
  double lmax = tot_load; //
  // fixme: do all reductions in one go
  mpi_all_reduce (tot_load, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (nmax, MPI_LONG, MPI_MAX);
  mpi_all_reduce (nmin, MPI_LONG, MPI_MIN);
  mpi_all_reduce (lmax, MPI_DOUBLE, MPI_MAX); // max load
  mpi_all_reduce (grid->tn, MPI_LONG, MPI_SUM);
  mpi_all_reduce (grid->maxdepth, MPI_INT, MPI_MAX);
  if (mpi.leaves)
    nt = grid->tn;
  else
    mpi_all_reduce (nt, MPI_LONG, MPI_SUM);

  assert (tot_load == nt); // todo: debug assertion, passed if w is unity
  double target_load = tot_load/npe(); // ideal load on each processor

  long ne = max(1, nt/npe());

  // todo: minimum number of cells?
  if (ne < mpi.min) {
    mpi.npe = max(1, nt/mpi.min);
    ne = max(1, nt/mpi.npe);
  }
  else
    mpi.npe = npe();

  // todo: early stopping? 99% efficency? 1% load? Maybe is never reached
  if (nmax - nmin <= 1)
    return false;

  if (pid() == 0) 
    fprintf (stderr, "efficency = %g\n", tot_load/npe()/lmax);
  if ((tot_load/npe())/lmax > 0.99)
    return false;

  scalar newpid[];
  double zn = z_indexing (newpid, mpi.leaves);
  if (pid() == 0)
    assert (zn + 1 == nt);

  FILE * fp = NULL;
#ifdef DEBUGCOND
  extern double t;
  if (DEBUGCOND)
    fp = lfopen ("bal", "w");
#elif DEBUG_MPI
  fp = lfopen ("bal", "w");
#endif

  // compute new pid, stored in newpid[]
  bool next = false, prev = false;
  double running_load = 0.;
  int current_rank = 0;
  foreach_cell_all() {
    if (is_local(cell)) {

      //int pid = current_rank;
      running_load += w[];
      if (running_load >= target_load && current_rank < npe() - 1) {
        current_rank += 1;
        running_load = 0.;
      }
      fprintf (stderr, "id: %g, rl: %g, pid: %d", newpid[], running_load, current_rank); 

      int pid = balanced_pid (newpid[], nt, mpi.npe);
      pid = clamp (pid, cell.pid - 1, cell.pid + 1);
      if (pid == pid() + 1)
        next = true;
      else if (pid == pid() - 1)
        prev = true;
      NEWPID()->pid = pid + 1;
      NEWPID()->leaf = is_leaf(cell);
      NEWPID()->prolongation = is_prolongation(cell);
      if (fp)
        fprintf (fp, "%g %g %d %d newpid\n", x, y, NEWPID()->pid - 1, cell.pid);
    }
    else
      newpid[] = 0;
  }
  for (int l = 0; l <= depth(); l++)
    boundary_iterate (level, {newpid}, l);

#ifdef DEBUGCOND
  extern double t;
  NOT_UNUSED(t);
  if (DEBUGCOND) {
    char name[80];
    sprintf (name, "colls-before-%d", pid());
    FILE * fp = fopen (name, "w");
    output_cells (fp);
    fclose (fp);

    sprintf (name, "pid-before-%d", pid());
    fp = fopen (name, "w");
    foreach_cell() {
      fprintf (fp, "%g %g %g %d %d %d %d\n",
          x, y, z, cell.pid, NEWPID()->pid - 1,
          NEWPID()->leaf, is_leaf(cell));
      if (NEWPID()->leaf)
        assert (is_leaf(cell));
    }
    fclose (fp);  
  }
#endif // DEBUGCOND

  Array * anext = next ? neighborhood (newpid, pid() + 1, fp) : array_new();
  Array * aprev = prev ? neighborhood (newpid, pid() - 1, fp) : array_new();

  if (fp)
    fflush (fp);

  check_flags();

  // send mesh to previous/next process
  MPI_Request rprev[2], rnext[2];
  if (pid() > 0)
    send_tree (aprev, pid() - 1, rprev);
  if (pid() < npe() - 1)
    send_tree (anext, pid() + 1, rnext);

  // receive mesh from next/previous process
  if (pid() < npe() - 1)
    receive_tree (pid() + 1, newpid, fp);
  if (pid() > 0)
    receive_tree (pid() - 1, newpid, fp);

  /* check that mesh was received OK and free send buffers */
  if (pid() > 0)
    wait_tree (aprev, rprev);
  array_free (aprev);
  if (pid() < npe() - 1)
    wait_tree (anext, rnext);
  array_free (anext);

  if (fp)
    fflush (fp);

  // set new pids
  int pid_changed = false;
  foreach_cell_all() {
    if (cell.pid >= 0) {
      if (is_newpid()) {
        if (fp)
          fprintf (fp, "%g %g %g %d %d %d %d %d new\n",
              x, y, z, NEWPID()->pid - 1, cell.pid,
              is_leaf(cell), cell.neighbors, NEWPID()->leaf);
        if (cell.pid != NEWPID()->pid - 1) {
          cell.pid = NEWPID()->pid - 1;
          cell.flags &= ~(active|border);
          if (is_local(cell))
            cell.flags |= active;
          pid_changed = true;
        }
        if (NEWPID()->leaf && !is_leaf(cell) && cell.neighbors)
          coarsen_cell_recursive (point, NULL);
      }
      else if (level > 0 && ((NewPid *)&coarse(newpid))->leaf)
        cell.pid = aparent(0).pid;
    }
    // cleanup unused prolongations
    if (!cell.neighbors && allocated_child(0)) {
      if (fp)
        fprintf (fp, "%g %g %g %d %d freechildren\n",
            x, y, z, NEWPID()->pid - 1, cell.pid);
      free_children (point);
    }
  }

  if (tree->dirty || pid_changed) {
#if 1
    // update active cells: fixme: can this be done above
    foreach_cell_post (!is_leaf (cell))
      if (!is_leaf(cell) && !is_local(cell)) {
        unsigned short flags = cell.flags & ~active;
        foreach_child()
          if (is_active(cell)) {
            flags |= active; break;
          }
        cell.flags = flags;
      }
#endif
    flag_border_cells(); // fixme: can this be done above?
    pid_changed = true;
  }

  if (fp)
    fclose (fp);

  mpi_all_reduce (pid_changed, MPI_INT, MPI_MAX);
  if (pid_changed)
    mpi_boundary_update_buffers();

  return pid_changed;
}


