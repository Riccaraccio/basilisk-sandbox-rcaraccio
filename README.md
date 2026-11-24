# Riccardo Caraccio
I am a PhD student in the [CRECK modeling group](https://www.creckmodeling.polimi.it/) at [Politecnico di Milano](https://www.polimi.it/).
My research is the multiscale modeling of biomass pyrolysis and combustion.
This sandbox holds my code for the flow in porous media and for the conversion of biomass particles.

Many thanks to [Edoardo Cipriano](https://basilisk.fr/sandbox/ecipriano/README).
His work gave this sandbox much of its ideas and its general structure.

The phase-change cases (biomass pyrolysis and combustion) need the [OpenSMOKE++ framework](https://www.opensmoke.polimi.it/).
This library gives the reaction rates, the physical properties of the species, and an efficient time integration of the chemistry.
OpenSMOKE++ is a C++ library, so Edoardo Cipriano and I wrote [a C interface to OpenSMOKE++](https://github.com/edocipriano/OpenSMOKEppInterface).

For more information, questions, or errors, contact me at "riccardo.caraccio@polimi.it".

## Publications
~~~bib
@article{CARACCIO2026174297,
  title = {Multiphase modeling of anisotropic biomass particle pyrolysis accounting for particle deformation and coupled gas-phase dynamics},
  journal = {Chemical Engineering Journal},
  volume = {531},
  pages = {174297},
  year = {2026},
  issn = {1385-8947},
  doi = {10.1016/j.cej.2026.174297},
  url = {https://www.sciencedirect.com/science/article/pii/S1385894726017560},
  author = {Riccardo Caraccio and Edoardo Cipriano and Alessio Frassoldati and Tiziano Faravelli},
  keywords = {Biomass, Volume of fluid, Porous, Shrinking, Pyrolysis, CFD}
}
~~~

## Requirements
* [Basilisk](http://basilisk.fr). Set `$BASILISK` to `basilisk/src`.
* [OpenSMOKE++](https://www.opensmoke.polimi.it/) and its [C interface](https://github.com/edocipriano/OpenSMOKEppInterface), for the phase-change cases only.
* [Vofi](https://github.com/VOFTracking/Vofi), for the geometric initialization of the VOF field.
* [GSL](https://www.gnu.org/software/gsl/), for the nonlinear solvers of the interface problems.

`run/Makefile` holds the include paths of these libraries.
If a library is in a different location on your machine, change `CFLAGS` in that file.

## Build and run
The cases are in `run/`, and the tests are in `test/`.
Build and run a case from its directory with the Basilisk test harness:

~~~bash
cd run
make porous-plug.tst
~~~

For a parallel run, use `CC='mpicc -D_MPI=4' make <case>.tst`.
The harness makes the directory `run/<case>/`, which holds the log, the dumps, and the plots.
Use `make <case>/plots` to make the figures that the source file declares.

The cases with OpenSMOKE++ need its library at run time.
Before you run such a case, set `LD_LIBRARY_PATH` to the `build/` directory of the C interface.
`test/Makefile` does this for the tests.

## Code structure
The headers in [src/](src/) add these models to the Basilisk solvers:

* [navier-stokes/centered-phasechange.h](src/navier-stokes/centered-phasechange.h): the centered Navier-Stokes solver with a gas mass source in the projection step. The source gives the volume expansion from the heterogeneous phase change.
* [darcy.h](src/darcy.h): the Darcy-Forchheimer penalization of the porous medium. The solver applies it implicitly, so it stays stable for very small permeabilities.
* [shrinking.h](src/shrinking.h): the shrinkage of the solid particle and the partition of the solid conversion between shrinkage and porosity.
* [multicomponent-varprop.h](src/multicomponent-varprop.h): the transport of the species and of the temperature in the gas phase and in the solid phase, with variable properties.
* [chemistry.h](src/chemistry.h) and [reactors.h](src/reactors.h): the gas-phase and solid-phase kinetics, with the OpenSMOKE++ stiff solver.
* [radiation.h](src/radiation.h): the radiation from the surface of the particle.
* [balances.h](src/balances.h): the mass and energy balances, for diagnostics.

### Acceleration of the gas chemistry
The stiff solve of the gas chemistry is the most expensive part of a combustion case.
These options decrease its cost:

* A frozen-cell gate skips the stiff solve of a gas cell that does not react.
* [binning.h](src/binning.h) groups cells with a similar state into bins, and the solver integrates one reactor for each bin ([Goldin et al. 2009](https://doi.org/10.1080/13647830903154542)).
* [drg.h](src/drg.h) reduces the mechanism with the Directed Relation Graph at the start of each stiff solve. This gives a speed-up only with large mechanisms.

## Simulations
This section gives a short description of the main simulations in this sandbox.

### Porous media
These cases test the [Darcy-Forchheimer penalization term](src/darcy.h).
This term adds the effect of a porous medium to the one-field Navier-Stokes equations.

* [Wake behind a porous cylinder](run/porous-cylinder.c)
* [Flow through a porous plug](run/porous-plug.c)
* [Flow in a partially porous channel](run/porous-channel.c)
* [Poiseuille flow through a porous section of a channel](run/porous-tube.c)

### Biomass pyrolysis
* [Mass conservation with phase change in a porous medium](run/massb.c): a comparison with an analytical solution for different grids and porosity models.
* [Pyrolysis of a biomass sphere in a temperature profile](run/corbetta.c): the temperatures in the particle and the gas species yields, with the data of Corbetta et al.
* [Pyrolysis of a biomass particle in a nitrogen flow](run/particle-flow.c): the mass loss and the shrinkage from 673 K to 973 K, with the data of Huang et al.
* [Pyrolysis of a biomass cylinder in a nitrogen flow](run/cylinder-flow.c): the mass loss and the radial and axial shrinkage at 723 K.

### Biomass combustion and gasification
* [Combustion of a biomass pellet in hot air (Fatehi)](run/fatehi-combustion.c): the solid mass, and the temperature and the H2O along the path of a laser beam.
* [Combustion of a biomass sphere in hot air (Lu)](run/lu-combustion.c): the solid mass, and the temperature at the core and at the surface.
* [Combustion of small biomass particles with different shapes (Glarborg)](run/glarborg-particle.c): one sphere and four cylinders with the same volume.
* [Gasification of a biomass cylinder in CO2 and O2](run/gasification.c): the solid mass at 1408 K.

### Gas-phase flames
* [Axisymmetric H2/N2-air coflow flame (Toro)](run/laminarflame.c): a validation of the gas-phase solver and of the chemistry acceleration, without a solid phase.

### Marangoni flow in a drop
* [Effective diffusivity from a Marangoni recirculation](run/marangoni-deff.c): the increase of the uptake of a tracer in a drop.
* [Sweep in the Marangoni and Schmidt numbers](run/marangoni-sweep.c), and [one point of the sweep for a cluster job array](run/marangoni-point.c).

## Tests
The tests in [test/](test/) check one part of the code each, mostly without flow.
They check the gas and solid chemistry, the binning, the DRG reduction, the frozen-cell gate, the Strang split, the gas expansion source, and the conservation of some interface operations.
Run a test from its directory:

~~~bash
cd test
make chemistry.tst
~~~

If the file `test/<test>.ref` exists, the harness compares the log of the test with it.

## Editor support for Basilisk C
I also maintain [basilisk-editor](https://github.com/Riccaraccio/basilisk-editor), which provides editor support for Basilisk C in vim/neovim, Emacs, VS Code and nano, with configurations for Helix, Sublime Text and Kate.
It recognises Basilisk files and highlights Basilisk constructs (`foreach` loops, events, fields, literate comments).
It also shows qcc errors as you type, and offers completion, documentation on hover and go to definition into the Basilisk headers.
It is still in beta (mostly developed with the help of an AI assistant), so feedback and bug reports are very welcome on the [GitHub issues page](https://github.com/Riccaraccio/basilisk-editor/issues).

## Co-authored publications
This is a short list of publications on related topics.
My other research is on my [Google Scholar page](https://scholar.google.com/citations?user=haLeMlYAAAAJ&hl=en&oi=sra).

~~~bib
@hal{cipriano2025, hal-04962085}
~~~
