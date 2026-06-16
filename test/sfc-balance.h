/**
## new partition policy

given a leaf's z-order index, the global leaf count, the
number of processes and the cumulative-weight array, return the
process id
* *cw*: cumulative-weights up to this cell.
* *tw*: total weight of the overall process.
* *nproc*: number of available processors.
*/

static int new_balanced_pid (double cw, double tw, int nproc) {
  double il = tw/nproc; // ideal load per process
  int idx = il > 0. ? (int) (cw/il) : 0;
  return min (idx, nproc - 1);
}

/**
## cumulative weight

`cf[i]` is the cumulative sum of `w` over leaves, in the same
order they are visited by the partitioning loop. */

void compute_cf (double* cf, long nt, (const) scalar w) {
  long idx = 0;
  double acc = 0.;
  foreach_leaf() {
    acc += w[];
    cf[idx++] = acc;
  }
}

/**
## per-process cell count and load 
performance monitor function that returns the number of cells 
and the load in each processor*/

void balance_score (long* counter, double* load, (const) scalar w) {
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
## static partitioning using the new policy

copy of `mpi_partitioning()` with the new partitioning method*/

trace
void new_mpi_partitioning ((const) scalar w = unity) {
  prof_start ("new_mpi_partitioning");

  long nt = 0;
  foreach (serial) nt++;
  
  double tl = 0.;
  foreach (serial)
    tl += w[];

  double cw = 0.;
  tree->dirty = true;
  foreach_cell_post (is_active (cell))
    if (is_active (cell)) {
      if (is_leaf (cell)) {

        cell.pid = new_balanced_pid (cw, tl, npe());
        cw += w[];

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

/**
# *z_weights()*: fills *cw* with the cumulative z-ordering weights.
   
if `leaves` is `true` only leaves are used, otherwise all active
cells. 

on the master process (`pid() == 0`), the function returns the
total load (and -1 on all other processes).

on a single processor, we would just need something
like (for leaves), given a weight scalar field w;

~~~literatec
scalar w[];
double rw = 0; // running weight
foreach() {
  cw[] = rw;   // weight of all leaves *before* this one
  rw += w[];
}
~~~

in parallel, this is a bit more difficult. */

trace
double z_weights (scalar cw, (const) scalar w, bool leaves)
{
  /**
  This function fills 'cw' with the cumulative-weight computed along
  the z-order for each cell given a weight scalar field 'w' (if 'leaves'
  is true, only leaves are counted).

  we first compute the weight of each subtree, stored in 'sw'. this is
  the weight analog of `size` in z_indexing(): we keep 'w' intact (the
  per-cell own weight) and accumulate the subtree sums into 'sw' via a
  (parallel) restriction from fine to coarse. a separate field is
  required because the indexing loop below needs both the parent's own
  weight and each child's subtree weight at the same time. */

  scalar sw[];                             // is this field needed?
  foreach()
    sw[] = w[];                            // own weight (set on leaves)

  boundary_iterate (restriction, {sw}, depth());
  for (int l = depth() - 1; l >= 0; l--) {
    foreach_coarse_level(l) {
      double sum = leaves ? 0. : w[];      // parent's own contribution
      foreach_child()
        sum += sw[];
      sw[] = sum;
    }
    boundary_iterate (restriction, {sw}, l);
  }

  /**
  The total load is the weight of the entire tree (i.e. the value of
  `sw` in the root cell on the master process). */

  double tl = -1.;
  if (pid() == 0)
    foreach_level(0, serial)
      tl = sw[];

  /**
  We now push the cumulative offset down the tree, exactly as
  z_indexing() does with the integer index. 'cw[]' on a cell is the
  cumulative weight of its subtree; the first child starts at the 
  parent's offset (plus the parent's own weight when indexing all cells), 
  and each subsequent sibling is offset by the preceding sibling's 
  subtree weight 'sw[]'. */

  // seed the root: master cell starts at 0
  foreach_level(0)
    cw[] = 0.;

  for (int l = 0; l < depth(); l++) {
    boundary_iterate (restriction, {cw}, l); // sync cw across mpi processes
    foreach_cell() {
      if (level == l) {
        if (is_leaf(cell)) {
          if (is_local(cell) && cell.neighbors) {
            double i = cw[]; // ghost children share the offset
            foreach_child()
              cw[] = i;
          }
        }
        else { // not leaf
          bool loc = is_local(cell);
          if (!loc)
            foreach_child()
              if (is_local(cell)) {
                loc = true; break;
              }
          if (loc) {
            double i = cw[] + (leaves ? 0. : w[]); // parent's own slot
            foreach_child() {
              cw[] = i;
              i += sw[]; // child's subtree weight
            }
          }
        }
        continue; // level == l
      }
      if (is_leaf(cell))
        continue;
    }
  }
  boundary_iterate (restriction, {cw}, depth());

  return tl;
}

trace
bool new_balance((const) scalar w = unity) {

  if (npe() == 1)
    return false;

  assert (sizeof(NewPid) == sizeof(double));

  check_flags();

  long nl = 0, nt = 0;
  foreach_cell() {
    if (is_local(cell)) {
      nt++;
      if (is_leaf(cell))
        nl++;
    }
    if (is_leaf(cell))
      continue;
  }

  grid->n = grid->tn = nl;
  grid->maxdepth = depth();
  long nmin = nl, nmax = nl;

  // fixme: do all reductions in one go
  mpi_all_reduce (nmax, MPI_LONG, MPI_MAX);
  mpi_all_reduce (nmin, MPI_LONG, MPI_MIN);
  mpi_all_reduce (grid->tn, MPI_LONG, MPI_SUM);
  mpi_all_reduce (grid->maxdepth, MPI_INT, MPI_MAX);

  if (mpi.leaves)
    nt = grid->tn;
  else
    mpi_all_reduce (nt, MPI_LONG, MPI_SUM);

  long ne = max(1, nt/npe());

  if (ne < mpi.min) {
    mpi.npe = max(1, nt/mpi.min);
    ne = max(1, nt/mpi.npe);
  }
  else
    mpi.npe = npe();

  scalar newpid[];
  double tl = z_weights (newpid, w, mpi.leaves);
  // We need to know the total load on all processors
  mpi_all_reduce (tl, MPI_DOUBLE, MPI_MAX);

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
  foreach_cell_all() {
    if (is_local(cell)) {
      int pid = new_balanced_pid (newpid[], tl, mpi.npe);
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

/**
## Automatic Load Balancing weights with TRACE
If 'TRACE' is enabled we can use it to compute the weight field automatically.
*/

#ifdef LB_AUTO
# if TRACE < 2
#  error "LB_AUTO requires an MPI build with tracing: compile with -D_MPI=<n> -DTRACE=2"
# else

/**
Traced functions where a rank is genuinely blocked waiting for other ranks.
Names not present in Trace.index (e.g. on a multigrid where mpi_waitany is
never called) simply never match, so listing extras is harmless.
*/
static const char * mpi_wait_funcs[] = {
  "mpi_waitany",      // point-to-point halo receive waits
  "rcv_pid_wait",     // send-completion waits (MPI_Wait)
  "mpi_all_reduce0",  // collective all-reduce / global synchronization
  NULL
};

/**
Helper function to extract the cumulative time of a list of specific functions
*/
static double trace_self_time (const char ** names) {
  double s = 0.;
  int len = Trace.index.len/sizeof(TraceIndex);
  TraceIndex * t = (TraceIndex *) Trace.index.p;
  for (int i = 0; i < len; i++, t++)
    for (const char ** n = names; *n; n++)
      if (!strcmp (t->func, *n)) {
        s += t->self;
        break;
      }
  return s;
}

/**
Fills weights with the values computed trough the tracing algorithm.
Each core has a total time t_wall = t_busy + t_wait. We want to extract 
t_busy = t_wall - t_wait. Therfore, we compute the elapsed time of the 
wait functions between each trace_weight function call.
*/
void trace_weights (scalar weights) {
  static double t_prev = -1., wait_prev = 0.;
  struct timeval tv; gettimeofday (&tv, NULL); // Could use Trace.t0

  double now = tv.tv_sec + tv.tv_usec/1e6;
  if (t_prev < 0) { // At the first execution we just initialize t_prev and wait_prev
    t_prev = now;
    wait_prev = trace_self_time (mpi_wait_funcs);
    foreach()
      weights[] = 1.; // Return uniform weights
    return;
  }

  double dt_wall = now - t_prev; t_prev = now;
  double wait_now = trace_self_time (mpi_wait_funcs);

  if (wait_now < wait_prev) 
    wait_prev = 0.;   // trace_print() may have reset counter

  double dt_wait = wait_now - wait_prev; wait_prev = wait_now;
  double t_busy = max (dt_wall - dt_wait, 0.);

  // Not exact per cell but should converge to the optimal split
  double per_cell = grid->n > 0 ? t_busy/grid->n : 0.;

  foreach()
    weights[] = per_cell + 1e-30;
}
# endif // TRACE < 2

# ifndef LB_ITER
#  define LB_ITER 1 // refresh weights every LB_ITER steps
# endif

#endif  // LB_AUTO

