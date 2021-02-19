#!/bin/sh

# Install macOS command line tools (fails if already present)
xcode-select --install

# Install package manager for macOS/Linux
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"

# Install required dependencies
brew install cmake gcc pkg-config

# Install some frequently used optional dependencies
brew install open-mpi gmp superlu suite-sparse eigen metis
