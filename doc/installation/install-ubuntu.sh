#!/bin/sh

# Update package manager
apt update

# Install required dependencies
apt install -y git cmake pkg-config g++ gfortran

# Install some frequently used optional dependencies
apt install -y libsuitesparse-dev openmpi-bin libopenmpi-dev libparmetis-dev libarpack++2-dev libsuperlu-dev

# Run the setup-dune script without root privileges
sudo -u $SUDO_USER ./doc/installation/setup-dune.sh
