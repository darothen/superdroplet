Fortran implementation of simple stochastic collision/coalescence algorithm, for testing purposes.

# Installation

A `CMakeLists.txt` with build recipes is available; it highlights dependencies, Fortran configurations, and external libraries needed for linking. For convenience we use Fortran bindings to the [GNU Scientific Library](https://www.gnu.org/software/gsl/) (fgsl); you might need to install them yourself and configure the *include_directories* in the provided cmake configuration.

The included `quick` script will run cmake for you and link an executable in this directory.

# Configuration / Running

You can configure a few parameters at the top of `sce.f90` in order to change the initial droplet size distribution, number of superdroplets, and simulation timestepping. Once changed, you'll need to rebuild the executable.

To run the compiled simulation, just execute the `sce` file.

A sequence of text output files will be generated with the schema **{timestep}_output.txt**. These files are normal CSV files with no header or index column; the first column is the superdroplet radius cubed, the second column is the multiciplity of that superdroplet.
