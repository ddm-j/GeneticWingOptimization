# GeneticWingOptimization

## Bespoke Evolutionary Algorithm for Wing Design

I wrote this package in my early development days. At the time I was not aware that some "wheels" of programming had already been made. So I set forth to re-make the wheel and my own evolutionary algorithm in Python. I have since been made aware of tailor-made GE packages like `DEAP`.

This package is meant to evolve a completely parametric aircraft wing in 3D space. The fitness function is the aerodynamic ratio of the wing, as computed via lift/drag coefficients. This package is extremely computationally intensive, and was originally designed for deployment on an HPC cluster; it may take days/weeks/months of computation time on a standard PC, depending on configuration.

## Simulation Setup

The meat of the computation is done using OpenFOAM 5.0. This is an opensource PDE solver that is commonly used for fluid-thermal CFD. The python code completely controls the work of the solver: creating case files, monitoring convergence, and calculating results.

## Dependencies

I am currently working on getting a neat set of dependencies for this old project. Here are some of the system considerations:

Linux OS (OSX/Win never been tested)
OpenFOAM 5.0 (Not been tested with newer OpenFOAM versions)
