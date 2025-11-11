# Riccardo Caraccio
I am a PhD student working at the [CRECK modeling group](https://www.creckmodeling.polimi.it/), [Politecnico di Milano](https://www.polimi.it/).
My research activity mainly focuses on the multiscale modeling of biomass pyrolysis and combustion.
In this sandbox you will find my code related to the simulation of flow in porous media and conversion of biomass particle.

Much love and thanks to the work done by [Edoardo Cipriano](https://basilisk.fr/sandbox/ecipriano/README) from which this works takes lots of inspiration and overall structure.

Please note that the the code realted to phasechange simualation ([Biomass pyrolysis](biomasspyro) and [combustion](biomasscomb)) requries integration with the [OpenSMOKE++ framework](https://www.opensmoke.polimi.it/).
This library is essential since we need to be able to readily compute reaction rates, species phyiscal properties and efficent chemistry time integration.
Additionally, given that the library is written in C++, [a C interface to OpenSMOKE++](https://github.com/edocipriano/OpenSMOKEppInterface) was developed by [Edoardo Cipriano](https://basilisk.fr/sandbox/ecipriano/README) and myself.

For further information, curiosities, or errors don't hesitate to contact me at "riccardo.caraccio@polimi.it"

## Publications
~~~bib
@article{caraccio2025,
  title={A volume-of-fluid model for biomass particle pyrolysis},
  author={Caraccio, Riccardo and Cipriano, Edoardo and Frassoldati, Alessio and Faravelli, Tiziano},
  journal={arXiv preprint arXiv:2510.17588},
  year={2025},
  month={November},
  http={https://arxiv.org/abs/2510.17588},
  DOI={10.48550/arXiv.2510.17588},
  pdf={https://arxiv.org/pdf/2510.17588.pdf}
}
~~~

## Simulations
A brief description of the main simulations reported in this sandbox.

### Porous media
Various simulations testing the implementation of the [Darcy-Forchheimer penalization term](src/darcy.h) that allows accounting for the presence of a porous medium in the one-field Navier-Stokes equations.

* [Wake behind a porous cylinder](run/porous-cylinder.c)
* [Flow obstructed by a porous plug](run/porous-plug.c)
* [Flow in a partially porous channel](run/porous-channel.c)

### Biomass pyrolysis [biomasspyro]
Coming soon!

### Biomass combustion [biomasscomb]
Coming soon!

## List of Co-Authored pubblications
Brief list of related topic on which I collaborated.
More research in different area can be found [here](https://scholar.google.com/citations?user=haLeMlYAAAAJ&hl=en&oi=sra)

~~~bib
@hal{cipriano2025, hal-04962085}
~~~
