/*

** TO BE REVISED **
OpenSMOKE rate and cell solution


*/

#include "osreactors.h"

event init (i = 0) {
  OpenSMOKE_InitODESolver ();
}

event cleanup (t = end) {
  OpenSMOKE_CleanODESolver ();
}

double GAMMA = 0.9;

#ifndef TURN_OFF_REACTIONS
event reaction (i++) {

	#ifdef SOLVE_TEMPERATURE
		odefunction batch = &batch_nonisothermal_constantpressure;
		unsigned int NEQ = NGS  + NSS + 1 + 1;
	#else
		odefunction batch = &batch_isothermal_constantpressure;
		unsigned int NEQ = NGS + NSS + 1 ;
	#endif

	double radius = sqrt (statsf(f).sum/pi);

	foreach(){
		if(f[]>(F_ERR)){

			////////////////////////////////////////////////////////
			//unpack the gas and solid fractions
			// foreach_elem(YGList, jj){
			// 	scalar YG = YGList[jj];
			// 	YG[] = (f[]<F_ERR) ? 0. : YG[]/f[];
			// }

			for(int jj=0; jj<NSS; jj++){
				scalar YS = YSList[jj];
				YS[] /= f[];
			}
      
			epsilon[] /= f[]; // check is not needed since f[] is already checked

			////////////////////////////////////////////////////////

			omega = 0.;
			double y0ode[NEQ]; 
			UserDataODE data;
			data.P = Pref + p[]; 
			data.T = T[];	
			data.rho = rho2;
			data.cp = cp2; 

			double gasmass[NGS];
			for (int jj=0; jj<NGS; jj++){
				scalar YG = YGList[jj];
				gasmass[jj] = YG[]*rho2*((1-epsilon[])*f[] + 1-f[]);
				y0ode[jj] = gasmass[jj];
			}

			double solidmass[NSS];
			for (int jj=0; jj<NSS; jj++){
				scalar YS = YSList[jj];
				solidmass[jj] = YS[]*rho1*epsilon[]*f[];
				y0ode[jj+NGS] = solidmass[jj];
			}

			if (sqrt(sq(x - 0.5*L0) + sq(y - 0.5*L0)) > radius*GAMMA){
				data.isinterface = 1;
			} else {
				data.isinterface = 0;
			}
			data.isinterface = 0;
			// if (f[] < (1-F_ERR)){
			// 	// y0ode[NGS+NSS] = f[];
			// 	// data.epsi = epsilon[];
			// 	data.isinterface = 1;
			// } else {
			// 	// y0ode[NGS+NSS] = epsilon[];
			// 	data.isinterface = 0; 
			// }

			y0ode[NGS+NSS] = epsilon[];
	
			#ifdef SOLVE_TEMPERATURE
				y0ode[NGS+NSS+1] = T[];
			#endif

			OpenSMOKE_ODESolver (batch, NEQ, dt, y0ode, &data);

			double totgasmass = 0;
			for (int jj=0; jj<NGS; jj++){
				totgasmass += y0ode[jj];
			}
			for (int jj=0; jj<NGS; jj++){
				scalar YG = YGList[jj];
				YG[] = y0ode[jj]/totgasmass;
			}

			double totsolidmass = 0;
			for (int jj=0; jj<NSS; jj++){
				totsolidmass += y0ode[jj+NGS];
			}
 
			for (int jj=0; jj<NSS; jj++){
				scalar YS = YSList[jj];
				YS[] = (totsolidmass<T_ERR) ? 0. : y0ode[jj+NGS]/totsolidmass;
			}

			// if (data.isinterface == 1) {
			// 	f[] = y0ode[NGS+NSS];
			// } else {
			// 	epsilon[] = y0ode[NGS+NSS];
			// }

      epsilon[] = y0ode[NGS+NSS];

			#ifdef SOLVE_TEMPERATURE
				T[] = y0ode[NGS+NSS+1];
			#endif

			// set MevapVal
			  // if (f[] < (1-F_ERR)){
			  //   coord m = mycs (point, f); // costante di piano
			  //   double alpha = plane_alpha (f[], m); // intercetta
			  //   coord prel; // centroide del piano, relativo al centro della cella
			  //   double segment = plane_area_center (m, alpha, &prel); // dimensione della cella, normalizzata per Delta
			  //   double rapporto = (segment > 0.) ? Delta/segment : 0.;
			  //   //mEvap[] = omega*rapporto; //kg/m^2/s
        //   //fprintf(stderr, "Mevap: %g\n", mEvap[]);
			  // }
			//fprintf(stderr,"Omega: %g\n",omega);
      double rhopf = rho1*epsilon[] + rho2*(1-epsilon[]);
			//drhodt[] = -(1/rhopf)*omega;
      ////////////////////////////////////////////////////////

      for(int jj=0; jj<NSS; jj++){
      	scalar YS = YSList[jj];
      	YS[] = YS[]*f[];
      }

      epsilon[] *= f[];
    }
	}

	boundary(YGList);
	boundary(YSList);

}
#endif


event graphs (t+=1) {
  double intsolidmass = 0.;
  extern double R0;

  foreach (reduction(+:intsolidmass)) {
    for (int jj=0; jj<NSS; jj++) {
      scalar YS = YSList[jj];
      intsolidmass += rho1*f[]*epsilon[]*YS[];
    }
  }

  double radius = sqrt (statsf(f).sum/pi);
  fprintf(stderr, "%g %g %g\n", t, intsolidmass/mass0, radius/R0); 
}