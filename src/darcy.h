/**
 * # Darcy and Forchheimer terms
 * This file implements the Darcy and Forchheimer terms for flow in porous media.
 * They allow to account for the resistance to flow due to the presence of a porous matrix.
 * The implementation is based on the Ergun equation, which combines both Darcy's law
 * and Forchheimer's correction for high flow rates.
*/

#ifndef F_ERR
  #define F_ERR 1e-10
#endif

/**
 * Extern variables defined elsewhere
 * - porosity: scalar field representing the porosity of the medium
 * - f: scalar field representing the volume fraction of the presudo-phase
 * - rhoGv_S, muGv_S: scalar fields for variable density and viscosity (if VARPROP is defined)
 * - rhoG, muG: double values for constant density and viscosity (if VARPROP is not defined)
 * 
 *  In this step, we want to solve the following equation:
 *  $$
 *  \frac{1}{\rho_g}\frac{\partial \mathbf{u}}{\partial t} = 
 *  - \left[ \frac{\mu_g \epsilon_g \mathbf{v}_g}{\mathbf{Da}} + 
 *  \rho_g\frac{1.75}{\sqrt{150\epsilon_g^3}} \frac{\epsilon_g
 *  |\mathbf{v}_g|\mathbf{v}_g}{\sqrt{\mathbf{Da}}} \right]
 *  $$
 *  
 *  Note that the permeability tensor Da can reach very small values.
 *  Therfore, an implicit treatment of the Darcy and Forchheimer 
 *  terms must and is here implemented.
 */

extern scalar porosity;
extern scalar f;

#ifdef VARPROP
extern scalar rhoGv_S, muGv_S;
#else
extern double rhoG, muG;
#endif

/**
 * Da permeability tensor
 * Here defined as a constant very small value to represent a highly resistive 
 * medium. Users can modify this to represent different porous media.
 * Units: m^2
 * The implementation allows for anisotropic permeability, different in each 
 * coordinate direction.
 */

coord Da = {1e-10, 1e-10};

/**
 * Viscous term event
 * After the viscous term is computed, this event modifies the velocity field
 * to account for the Darcy and Forchheimer resistance. 
 */

event viscous_term (i++) {
  foreach() {
    if (f[] > F_ERR) {
      double e = porosity[]/f[];
      double F = 1.75/sqrt (150*pow (e, 3));
      double Umag = sqrt (sq(u.x[]) + sq(u.y[]));

      double muGh, rhoGh;
      #ifdef VARPROP
      muGh = muGv_S[];
      rhoGh = rhoGv_S[];
      #else
      muGh = muG;
      rhoGh = rhoG;
      #endif

      foreach_dimension() {
        double A = (1./rhoGh)*(muGh*e/Da.x);   // Darcy term
        double B = F*e/sqrt(Da.x)*Umag;        // Forchheimer term
        u.x[] *= exp(-(A + B)*dt*f[]);
      }
    }
  }
}

/**
 * Previous implementation (commented out)
 * Here we used the default acceleration field 'a' to add the Darcy and Forchheimer
 * contributions explicitly. This appoach is more integrated with the existing framework but
 * less efficient due to the explicit nature of the terms.
 */

// Old approach, works but explicit
// event defaults (i = 0) {
//   if (is_constant(a.x)) {
//     a = new face vector;
//     foreach_face() {
//       a.x[] = 0.;
//       dimensional (a.x[] == Delta/sq(DT));
//     }
//   }
// }

// event acceleration (i++){
//   face vector av = a;
//   foreach_face() {
//     double ff = face_value(f, 0);
//     if (ff > F_ERR) {
//       double ef = face_value(porosity, 0);
//       double F  = 1.75/pow (150*pow (ef, 3), 0.5);

//       // Darcy contribution, weighted by the face fraction of the interface
//       av.x[] -= alpha.x[]/(fm.x[] + SEPS)* (muG*ef/Da.x) *uf.x[] *ff; 

//       // Forcheimer contribution
//       av.x[] -= alpha.x[]/(fm.x[] + SEPS)* (F*ef*rhoG/pow(Da.x,0.5)) *fabs(uf.x[])*uf.x[] *ff;
//     }
//   }
// }