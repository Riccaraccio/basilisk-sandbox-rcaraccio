// We want to define a macro that returns the levelset function of a superquadric
// The superquadric is defined by the following equation:
// (x/a)^m + (y/b)^m = 1
// where m is the shape parameter, a is the radius in the x direction, and b is the radius in the y direction
// The levelset function is defined as the distance from the surface of the superquadric

macro double superquadric (double x, double y, 
                            double m, double a, double b,
                            double xc = 0., double yc = 0.)
{
    return 1. - pow((x-xc*L0)/a, m) - pow((y-yc*L0)/b, m);
}

// macro double rectangle (double x, double y, 
//                     double a, double b, 
//                     double xc = 0., double yc = 0.)
// {
//   double distancex = fabs(x - xc*L0) - a;
//   double distancey = fabs(y - yc*L0) - b;

//   double outside_distance = max(max(distancex, distancey), 0.); //FIXME
//   double inside_distance = min(max(distancex, distancey), 0.);
//   return outside_distance + inside_distance;
// }