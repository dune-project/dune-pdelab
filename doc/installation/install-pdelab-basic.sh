#!/bin/sh

# Update package manager
sudo apt update

# Install required dependencies
sudo apt install -y git cmake pkg-config g++ gfortran

# Install some frequently used optional dependencies
sudo apt install -y libsuitesparse-dev openmpi-bin libopenmpi-dev libparmetis-dev libarpack++2-dev libsuperlu-dev


# Create new directory to build in
mkdir my-dune
cd my-dune

# Get DUNE modules via git
git clone https://gitlab.dune-project.org/core/dune-common.git
git clone https://gitlab.dune-project.org/core/dune-geometry.git
git clone https://gitlab.dune-project.org/staging/dune-uggrid.git
git clone https://gitlab.dune-project.org/core/dune-grid.git
git clone https://gitlab.dune-project.org/core/dune-localfunctions.git
git clone https://gitlab.dune-project.org/staging/dune-typetree.git
git clone https://gitlab.dune-project.org/core/dune-istl.git
git clone https://gitlab.dune-project.org/staging/dune-functions.git
git clone https://gitlab.dune-project.org/pdelab/dune-pdelab.git

# Choose a release version for all modules
./dune-common/bin/dunecontrol git checkout releases/2.7

# Build dune-pdelab and its dependencies
./dune-common/bin/dunecontrol --builddir=$PWD/build --module=dune-pdelab all


# Let's verify our build was successful by building and running a simple test
cd build/dune-pdelab/dune/pdelab/test
make testpoisson
./testpoisson
