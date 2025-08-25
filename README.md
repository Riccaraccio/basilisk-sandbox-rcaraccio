# Basilisk Sandbox for Gas-Solid Reacting Media Simulations

This repository contains simulation code for modeling gas-solid reacting media using Basilisk, a free software program for solving partial differential equations on adaptive Cartesian meshes.

## Overview

This collection of simulation codes implements advanced numerical models for studying complex physicochemical phenomena in gas-solid systems, including:

- **Multicomponent mass transport** with variable properties
- **Chemical reactions** using detailed kinetic mechanisms via OpenSMOKE++ integration
- **Heat transfer** with temperature-dependent properties
- **Porous media flow** using Darcy's law
- **Interface tracking** with shrinking particle models
- **Adaptive mesh refinement** for computational efficiency

The simulations are particularly suited for applications in:
- Biomass pyrolysis and gasification
- Solid fuel combustion
- Catalytic processes
- Packed bed reactors
- Particle devolatilization

## Features

### Core Capabilities
- **Multicomponent Gas Phase**: Support for multiple gas species with variable properties
- **Solid Phase Reactions**: Detailed solid-phase chemistry with species tracking
- **Interface Dynamics**: Moving boundaries with mass transfer across interfaces
- **Adaptive Grids**: Automatic mesh refinement based on solution gradients
- **Mass & Energy Balances**: Comprehensive conservation tracking and validation
- **Flexible Geometries**: Support for axisymmetric, cylindrical, and slab configurations

### Physical Models
- **Chemistry Integration**: OpenSMOKE++ for detailed chemical kinetics
- **Transport Phenomena**: Fick's law with corrections, variable diffusivity
- **Heat Transfer**: Conduction, convection, and interface heat exchange
- **Porous Flow**: Darcy velocity with porosity evolution
- **Radiation**: Interface radiation modeling (optional)

## Repository Structure

```
├── src/                          # Source code modules
│   ├── multicomponent-varprop.h  # Multicomponent transport with variable properties
│   ├── balances-interface.h      # Mass balance tracking and validation
│   ├── chemistry.h               # Chemical reaction integration
│   ├── darcy.h                   # Porous media flow
│   ├── shrinking.h               # Moving boundary conditions
│   ├── intgrad.h                 # Interface gradient calculations
│   ├── opensmoke.h               # OpenSMOKE++ integration
│   └── ...                       # Additional utility modules
├── run/                          # Simulation cases
│   ├── particle-flow.c           # Particle flow simulation
│   ├── ancacouce.c              # Specific case study simulation
│   ├── cylinder-flow.c          # Cylindrical geometry flow
│   ├── slab.c                   # Slab geometry simulation
│   └── ...                      # Additional test cases
├── data/                        # Experimental validation data
│   ├── corbetta/                # Corbetta et al. experimental data
│   ├── cylinder-flow/           # Cylinder flow validation data
│   └── ...                     # Additional datasets
└── README.md                    # This file
```

## Dependencies

### Required Software
- **Basilisk CFD**: The main computational framework ([basilisk.fr](http://basilisk.fr))
- **OpenSMOKE++**: Chemical kinetics library for detailed reactions
- **GSL**: GNU Scientific Library for numerical solvers
- **SUNDIALS**: Suite for nonlinear and differential/algebraic equation solvers

### Environment Variables
Ensure the following environment variables are set:
```bash
export OPENSMOKE_INTERFACE=/path/to/opensmoke/interface
export BASILISK=/path/to/basilisk
```

### System Requirements
- C compiler (GCC recommended)
- OpenMPI (for parallel simulations)
- GSL development libraries
- SUNDIALS development libraries

## Installation

1. **Install Basilisk**: Follow instructions at [basilisk.fr](http://basilisk.fr)

2. **Install OpenSMOKE++**: Obtain and compile OpenSMOKE++ with interface libraries

3. **Install system dependencies** (Ubuntu/Debian):
   ```bash
   sudo apt-get install libgsl-dev libsundials-dev
   ```

4. **Clone this repository**:
   ```bash
   git clone https://github.com/Riccaraccio/basilisk-sandbox-rcaraccio.git
   cd basilisk-sandbox-rcaraccio
   ```

## Usage

### Basic Compilation
Basilisk simulations are typically compiled and run in a single command:

```bash
# Compile and run a simulation
$BASILISK/qcc -O2 -Wall run/particle-flow.c -o particle-flow
./particle-flow
```

### Example Simulations

#### 1. Particle Flow Simulation
```bash
# Biomass particle flow with detailed chemistry
$BASILISK/qcc -O2 run/particle-flow.c -o particle-flow
./particle-flow
```

#### 2. Cylinder Flow Case
```bash
# Flow around cylindrical particles
$BASILISK/qcc -O2 run/cylinder-flow.c -o cylinder-flow
./cylinder-flow
```

#### 3. Slab Geometry
```bash
# One-dimensional slab configuration
$BASILISK/qcc -O2 run/slab.c -o slab
./slab
```

### Simulation Parameters
Most simulations can be customized by modifying parameters at the top of the source files:

- **Grid resolution**: `maxlevel` parameter
- **Physical properties**: Material properties (density, thermal conductivity, etc.)
- **Boundary conditions**: Inlet velocities, temperatures, compositions
- **Chemistry**: Kinetic mechanism folder path

### Output Files
Simulations typically generate:
- `balances-*`: Mass and energy balance files
- `OutputData-*`: Detailed field data
- Various field dumps for post-processing

## Key Simulation Cases

### Particle Flow (`particle-flow.c`)
- **Application**: Single particle devolatilization in crossflow
- **Features**: Multicomponent transport, variable properties, shrinking
- **Validation**: Experimental data comparison available

### Ancacouce Case (`ancacouce.c`)
- **Application**: Specific biomass pyrolysis study
- **Features**: Detailed chemistry, mass balances, adaptive grid
- **Geometry**: Superquadric particle shape

### Cylinder Flow (`cylinder-flow.c`)
- **Application**: Flow and heat transfer around cylindrical particles
- **Features**: Radiation, convection, mass transfer
- **Validation**: Comparison with experimental correlations

## Validation Data

The `data/` directory contains experimental datasets for validation:
- **Corbetta et al.**: Single particle pyrolysis experiments
- **Cylinder flow**: Heat and mass transfer correlations
- **Various datasets**: Additional validation cases

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-model`)
3. Commit your changes (`git commit -am 'Add new model'`)
4. Push to the branch (`git push origin feature/new-model`)
5. Create a Pull Request

### Code Style
- Follow existing naming conventions
- Document new modules with appropriate headers
- Include validation cases for new models
- Ensure mass/energy conservation in new implementations

## References

Key publications related to this work:
- Basilisk documentation: [basilisk.fr](http://basilisk.fr)
- OpenSMOKE++: Cuoci et al., Computer Physics Communications (2015)
- VOF methods: Popinet, Journal of Computational Physics (2009)
- Interface gradients: Fleckenstein & Bothe, Journal of Computational Physics (2015)

## License

This code is distributed under the same license as Basilisk (GPL). See individual files for specific licensing information.

## Support

For questions about the simulations or to report issues:
- Open an issue on GitHub
- Check the Basilisk documentation and forums
- Review the validation cases for examples

## Acknowledgments

This work builds upon the Basilisk CFD framework developed by Stéphane Popinet and the OpenSMOKE++ library developed by Alberto Cuoci and colleagues.
