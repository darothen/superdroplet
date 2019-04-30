C++ implementation of simple stochastic collision/coalescene algorithm, for testing purposes.

## Build Instructions

A standard CMake build script is supplied, as well as a script `quick` which automates the building process. In general, to build and run the model:

- Create an empty folder, **bld/**
- `cd` into **bld/** and execute the command:

  ``` shell
  cmake ..
  ```

  This will build an executable, `sce`, in place

Run the executable and you're ready to go.

## Model Configuration

All the necessary hooks to tweak the model are provided at the beginning of the *main()* method in **sce.cpp**. Note that only a single collision kernel (Golovin) has been implemented in this code so far.