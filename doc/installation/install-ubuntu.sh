#!/bin/sh

# We expect all commands to exit successfully
set -e

# Update package manager
apt update

# Install required dependencies
apt install -y git cmake pkg-config g++ gfortran

# Install some frequently used optional dependencies
apt install -y libsuitesparse-dev openmpi-bin libopenmpi-dev libparmetis-dev libarpack++2-dev libsuperlu-dev
