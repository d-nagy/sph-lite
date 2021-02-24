# sph-lite

sph-lite is Phase I of my master's level project with Tobias Weinzierl, "Mini-SPH".

The aim of Phase I is to develop a very simple test-bed SPH code that can be later integrated into the HPC software [Peano](http://www.peano-framework.org/index.php/peano-v-3/).

Usage (Linux only):
- Clone repository.
- Compile code with `make` in the repo folder (If you don't have clang++ then change CXX)
  - Alternatively run `make parallel` to compile with OpenMP parallelism
- Generate simulation initial conditions.
  - See `cases` dir for examples of python scripts which calculate the initical positions of particles.
  - Python script generates a `.case` file which is referenced by the setup file (see next bullet point)
- Create a setup file (see `sedov.dat`, `stationary.dat` and `dambreak.dat` for reference).
- Run simulation with `./sph-lite.exe <setup file>`


SPH Kernels implemented:
- Cubic spline

SPH schemes:
- Weakly Compressible SPH
- "Thermo SPH", using an ideal gas equation of state that incorporates internal energy

Particle types:
- Fluid
- Boundary - non-moving particles acting as rigid boundaries (see dambreak example)

Boundary conditions:
- Periodic
- Destructive (destroy particles which leave the domain)
