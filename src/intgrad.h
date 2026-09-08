/**
# Interface Gradients

The calculation of the interface gradients can be performed
using the method developed in [embed.h](/src/embed.h). This
module was extended from [gradients.h](/sandbox/ggennari/phase_change/gradients.h)
and it reports the same interface gradients calculation
implemented in [embed.h](/src/embed.h), together with the
vof-averaged interface gradients proposed by [Fleckenstein and Bothe](#fleckenstein2015volume):

![VOF-averaged normal gradient scheme](/src/figures/dirichlet_gradient.svg)

$$
\left(\dfrac{\partial f}{\partial \mathbf{n}_\Gamma}\right)
=
\left(c\dfrac{f_\Gamma - f_0}{d_0} + \left(1 - c\right)\dfrac{f_\Gamma - f_1}{d_1}\right)
$$
*/

/**
## Copy of embed

Part of the embedded gradient functions are copied here, with
different names in order to eventually use the gradients also
with the [embed.h](/src/embed.h) module.
*/

#define quadratic(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))

/**
## The two forms of the gradient, and why the second one exists

The stencil below reads the NEIGHBOURS of the cut cell along the normal. The
value of the cell itself is not in it. That is accurate, and it is also the
reason a sliver is not well posed: the flux that heats the cell does not
answer to the temperature of the cell, while the heat capacity and the face
conductance of the cell both go to zero with its phase fraction. The row of
such a cell has no diagonal and the solve does not bound it. Read the design
note at the head of `int-temperature-vofbc.h`.

So the function has a second form, the same one that `embed.h` uses in its
degenerate branch:

    grad = (bc - s[])/(d0*Delta)

This form reads the cell, so it puts a term on the diagonal and the row obeys
a maximum principle. It is first order, and it is the only form available when
the stencil is missing.

`force` selects between them.

- `force == false` — the historical behaviour, bit for bit. The accurate
  stencil, and a degenerate cell returns a zero gradient, which silently
  drops its flux.
- `force == true` — the first-order form, everywhere.

Neither form returns the two parts separately, and neither takes a pointer.
`ebmgrad_bc()` recovers the split from the first form alone, because it is
affine in `bc`:

    grad(bc) = (bc - s[])/(d0*Delta)
    grad(1) - grad(0) = 1/(d0*Delta)

so one difference gives both the diagonal coefficient and the source part,
with no output argument. That matters: the stencil pass of `qcc` must know
which neighbours a point function reads, and it cannot follow a pointer or a
variable flag through a call.

`d0` is the distance from the centroid of the phase that carries `s` to the
interface, in units of `Delta`. `ebmgrad` computes it. A value of zero asks
for the estimate of `embed.h`, `fabs(p.x/n.x)`, which measures from the cell
centre. For a sliver the cell centre is inside the OTHER phase, so that
estimate is far too large and the floor of `1e-3` does all the work. Pass the
centroid distance whenever you have it. */

foreach_dimension()
static inline double concentration_gradient_x (Point point, scalar s, scalar cs, face vector fs,
					   coord n, coord p, double bc,
					   bool third, double * coef,
					   double d0, bool force)
{
  foreach_dimension()
    n.x = - n.x;
  double d[2], v[2] = {nodata,nodata};
  bool defined = !force;
  foreach_dimension()
    if (defined && !fs.x[(n.x > 0.)])
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
#if dimension == 2
      if (fs.x[i + (i < 0),j] && fs.y[i,j] && fs.y[i,j+1] &&
	  cs[i,j-1] && cs[i,j] && cs[i,j+1])
	v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
#else // dimension == 3
      double z = p.z + d[l]*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = fs.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
	if (!fs.y[i,j,k+m] || !fs.y[i,j+1,k+m] ||
	    !fs.z[i,j+m,k] || !fs.z[i,j+m,k+1] ||
	    !cs[i,j+m,k-1] || !cs[i,j+m,k] || !cs[i,j+m,k+1])
	  defined = false;
      if (defined)
	// bi-quadratic interpolation
	v[l] =
	  quadratic (z,
		     quadratic (y1,
				(s[i,j-1,k-1]), (s[i,j,k-1]), (s[i,j+1,k-1])),
		     quadratic (y1,
				(s[i,j-1,k]),   (s[i,j,k]),   (s[i,j+1,k])),
		     quadratic (y1,
				(s[i,j-1,k+1]), (s[i,j,k+1]), (s[i,j+1,k+1])));
#endif // dimension == 3
      else
	break;
    }
  if (v[0] == nodata) {

    /**
    The stencil is not available, or the caller forced this branch. */

    if (!force) {

      /**
      The caller did not ask for the first-order form, so it has nowhere to
      put a diagonal term. Keep the historical behaviour: report no flux.
      This DROPS the interface flux of the cell, and nothing reports it.
      Count these cells before you rely on them. */

      *coef = 0.;
      return 0.;
    }

    double dd = max (1e-3, (d0 > 0.) ? d0 : fabs(p.x/n.x));
    *coef = - 1./(dd*Delta);
    return bc/(dd*Delta);
  }

  /**
  For non-degenerate cases, the gradient can be  obtained using
  second-order, third-order or vof-averaged estimates. None of the three
  reads `s[]`, so none of them puts a term on the diagonal. */

  *coef = 0.;

  if (v[1] != nodata && third) { // third-order gradient
    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
  }
  else if (v[1] != nodata) {
    return (cs[]*(bc - v[0])/d[0] + (1. - cs[])*(bc - v[1])/d[1])/Delta; //vof-avg
  }
  return (bc - v[0])/(d[0]*Delta); // second-order gradient
}

double concentration_gradient (Point point, scalar s, scalar cs, face vector fs,
				coord n, coord p, double bc, bool third,
				double * coef, double d0, bool force)
{
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return concentration_gradient_x (point, s, cs, fs, n, p, bc, third,
                                       coef, d0, force);
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return concentration_gradient_x (point, s, cs, fs, n, p, bc, third,
                                       coef, d0, force);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return concentration_gradient_y (point, s, cs, fs, n, p, bc, third,
                                     coef, d0, force);
  return concentration_gradient_z (point, s, cs, fs, n, p, bc, third,
                                   coef, d0, force);
#endif // dimension == 3
  return nodata;
}

/**
## The thickness of the phase, seen from the interface

The distance from the centroid of one phase of a cut cell to the interface,
along the normal, in units of `Delta`. This is the length scale over which
that phase can hold a gradient.

Caution: pass the BOX-normalised `m` and the `alpha` that belongs to it, that
is, before any `normalize(&m)`. `plane_alpha` returns the alpha of the box
normal, and normalising `m` breaks the pair.

The centroid of the second phase follows from the first and from the centre
of the cell, `cL*pL + cG*pG = 0`, so one call to `plane_center` serves both
sides. */

static inline double interface_phase_distance (Point point, scalar fL,
                                               coord m, double alpha,
                                               bool inverse)
{
  double cL = fL[], cG = 1. - cL;
  if (cL <= 0. || cG <= 0.)
    return 0.;

  coord pc;
  plane_center (m, alpha, cL, &pc);

  if (inverse)
    foreach_dimension()
      pc.x *= - cL/cG;

  double mp = 0., mm = 0.;
  foreach_dimension() {
    mp += m.x*pc.x;
    mm += sq(m.x);
  }
  return (mm > 0.) ? fabs (alpha - mp)/sqrt(mm) : 0.;
}

/**
## *ebmgrad()*: high-level interface for the calculation of the interface gradients:
* *tr*: scalar fields whose gradients must be computed
* *fL*: liquid phase volume fraction (*f*)
* *fG*: gas phase volume fraction (*1 - f*)
* *fsL*: face fraction for *fL*, computed using [fracface.h](fracface.h)
* *fsG*: face fraction for *fG*, computed using [fracface.h](fracface.h)
* *inverse*: true if tracer is in gas phase, false otherwise
* *trint*: interface value
* *success*: deprecated (fixme)

`ebmgrad_force()` takes one more argument and `ebmgrad()` is the wrapper that
passes the historical value of it. Two entry points, and not one function with
a default argument, because the stencil pass of `qcc` needs the flag to be a
constant at each call site.

* *force*: take the first-order form that reads this cell, even where the
  accurate stencil exists. A caller uses this to repair a row that would
  otherwise have no diagonal.
*/

double ebmgrad_force (Point point,
                  scalar tr,
                  scalar fL,
                  scalar fG,
                  face vector fsL,
                  face vector fsG,
                  bool inverse,
                  double trint,
                  bool* success,
                  bool force)
{
  coord m = interface_normal (point, fL);
#if dimension == 2
  coord p = {0.,0.};
#else //dimension ==3
  coord p = {0.,0.,0.};
#endif

  double alpha = plane_alpha (fL[], m);
  plane_area_center (m, alpha, &p);

  /**
  The distance from the centroid of the phase that carries `tr` to the
  interface, along the normal, in units of `Delta`. Only the first-order form
  uses it, so compute it only when a caller may take that form.

  Caution: this must run BEFORE `normalize(&m)`. `plane_alpha` returns the
  `alpha` of the BOX-normalised normal, and `normalize` makes `m` a Euclidean
  unit vector, after which `alpha` no longer belongs to `m`. */

  double d0 = force ? interface_phase_distance (point, fL, m, alpha,
                                                tr.inverse) : 0.;

  normalize (&m);

#if dimension == 2
  coord n = {0., 0.};
#else
  coord n = {0., 0., 0.};
#endif

  if (tr.inverse) {
    foreach_dimension()
      n.x = -m.x;
  }
  else {
    foreach_dimension()
      n.x = m.x;
  }

  bool third = false;
#ifdef INTGRAD_3rd
  third = true;
#endif

  /**
  The kernel splits its answer, so put it back together: this entry point
  returns the WHOLE gradient. `coef` is zero unless the first-order form was
  taken, so the reconstruction is exact and costs nothing in the usual case.
  `plicbc.h` calls the kernel directly when it wants the two parts apart. */

  double c = 0., dirgrad = 0.;
  if (tr.inverse)
    dirgrad = concentration_gradient (point, tr, fG, fsG, n, p, trint, third,
                                      &c, d0, force);
  else
    dirgrad = concentration_gradient (point, tr, fL, fsL, n, p, trint, third,
                                      &c, d0, force);

  return dirgrad + c*tr[];
}

/**
The historical entry point. It returns the whole gradient and it never takes
the first-order form, so every call site that predates the split is unchanged,
bit for bit. */

double ebmgrad (Point point,
                  scalar tr,
                  scalar fL,
                  scalar fG,
                  face vector fsL,
                  face vector fsG,
                  bool inverse,
                  double trint,
                  bool* success)
{
  return ebmgrad_force (point, tr, fL, fG, fsL, fsG, inverse, trint, success,
                        false);
}

/**
## *ebmgrad_bc()*: the gradient a linear solve can use

This is the entry point of `INT_TEMP_VOFBC`. It returns the gradient split
for a solve, and it decides which of the two forms the cell needs. The design
note is at the head of `int-temperature-vofbc.h`; the short version is that
the accurate stencil leaves the row of a sliver with no diagonal, so such a
cell must take the first-order form that reads itself.

The flux of the cell is

    lambdah*aov*(return + (*coef)*tr[])

* *lambdah*: the conductivity of this phase, weighted along the normal.
* *aov*: the interface area of the cell divided by its volume, with the
  metric. The caller builds it, because the axisymmetric form differs.
* *ddiag*: the diagonal the row would have without the interface term, that
  is `theta/dt` plus the face conductance. The cap uses it. The criterion
  does not: see the comment in the body.
* *eta*: take the first-order form when the thickness of the phase in this
  cell is below `eta` times the reach of the accurate stencil.
* *cmax*: cap the interface conductance at `cmax*ddiag`. Pass 0 for no cap.
  The cap scales the returned pair together, so it stays a flux of the form
  `C*(trint - tr[])` and the value the cell relaxes to does not move.
* *branch*: what the cell did. 0 the accurate stencil, 1 the criterion fired,
  2 the stencil was missing. Add 10 when the cap fired.
* *dratio*: the ratio the criterion tests, so a run can report what it saw.
  Pass `NULL` if you do not report it.

The slope `h` is exact as a difference of two calls, because `ebmgrad` is
affine in the interface value. Do not finite-difference it. */

double ebmgrad_bc (Point point,
                   scalar tr,
                   scalar fL,
                   scalar fG,
                   face vector fsL,
                   face vector fsG,
                   bool inverse,
                   double trint,
                   double lambdah,
                   double aov,
                   double ddiag,
                   double eta,
                   double cmax,
                   double * coef,
                   int * branch,
                   double * dratio)
{
  bool ok = false;

  /**
  The slope of the accurate form in the interface value. `ebmgrad` is affine
  in that value, so the slope is exact as a difference of two calls. Do not
  finite-difference it.

  The slope is also a length: the accurate estimate is
  `(bc - v0)/(d1*Delta)`, so `d1 = 1/(h*Delta)` is the distance to the point
  the stencil interpolates, in units of `Delta`. It is about one cell,
  whatever the phase fraction is.

  A cell with no stencil returns zero from both calls, so `h` is exactly
  zero. That is the test for the degenerate case, and it needs no separate
  output: no accurate form has a zero slope. */

  double h = ebmgrad (point, tr, fL, fG, fsL, fsG, inverse, 1., &ok)
           - ebmgrad (point, tr, fL, fG, fsL, fsG, inverse, 0., &ok);
  bool degenerate = (h == 0.);

  /**
  The criterion. `d0` is the thickness of this phase in the cell, and `d1` is
  the distance the stencil reaches to build the gradient. Their ratio says
  whether the interpolated value describes this phase at all.

  For a well resolved cut cell the two are comparable. For a sliver `d0` goes
  to zero while `d1` stays at about one cell, so the stencil reads a value
  from far outside the layer it claims to differentiate, and the flux it
  returns does not answer to the cell. Below `eta` the cell takes the
  first-order form instead, which reads itself.

  This is NOT the ratio of the interface conductance to the diagonal of the
  row. That ratio is bounded: for a convex cut the interface area is at most
  the sum of the face fractions of the phase, so the interface conductance is
  at most the face conductance, whatever the fraction. A sliver row is
  diagonally dominant and the solve converges; what it converges to is the
  wrong thing. Measured: the ratio sits at 0.33 from a gas fraction of 0.4
  down to 1e-6, while the answer overshoots the interface value by a factor
  of four at every one of them. See `test/intbc-sliver.c`. */

  coord m = interface_normal (point, fL);
  double d0 = interface_phase_distance (point, fL, m,
                                        plane_alpha (fL[], m), tr.inverse);
  double d1 = (h != 0.) ? 1./(fabs (h)*Delta) : 0.;
  double dr = (d1 > 0.) ? d0/d1 : 0.;
  if (dratio)
    *dratio = dr;

  bool force = degenerate || (dr < eta);

  if (!force) {
    *branch = 0;
    *coef = 0.;
    return ebmgrad (point, tr, fL, fG, fsL, fsG, inverse, trint, &ok);
  }

  /**
  The first-order form is `(bc - tr[])/(dd*Delta)`, so its slope in `bc` is
  `1/(dd*Delta)`. One difference therefore gives both parts of the split:

      coef = -1/(dd*Delta) = -(g(1) - g(0))
      G    =  bc/(dd*Delta) = trint*(g(1) - g(0))

  Two calls, both with a literal `force`, and no output pointer through which
  the stencil pass of `qcc` would have to follow the flag. */

  double gr = ebmgrad_force (point, tr, fL, fG, fsL, fsG, inverse, 1., &ok, true)
            - ebmgrad_force (point, tr, fL, fG, fsL, fsG, inverse, 0., &ok, true);

  double c = - gr, g = trint*gr;
  *branch = degenerate ? 2 : 1;

  /**
  The cap. Scale the source part and the diagonal part by the same factor, so
  the pair is still the flux `C*(trint - tr[])` with a smaller `C`, and the
  value the cell relaxes to does not move. */

  if (c != 0. && cmax > 0. && ddiag > 0.) {
    double cc = fabs (lambdah*c)*aov;
    if (cc > cmax*ddiag) {
      double sc = cmax*ddiag/cc;
      g *= sc;
      c *= sc;
      *branch += 10;
    }
  }

  *coef = c;
  return g;
}

/**
## References

~~~bib
@article{fleckenstein2015volume,
  title={A volume-of-fluid-based numerical method for multi-component mass transfer with local volume changes},
  author={Fleckenstein, Stefan and Bothe, Dieter},
  journal={Journal of Computational Physics},
  volume={301},
  pages={35--58},
  year={2015},
  publisher={Elsevier}
}
~~~
*/

