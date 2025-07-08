#!/bin/bash

which conda > /dev/null 2>&1 || { echo "conda not found..exiting"; exit 1;}

if [ $# -gt 0 ]; then
    QFORTE_CONDA_ENV=$1
else
    QFORTE_CONDA_ENV="qforte-default-env"
fi

echo "Initializing Conda...\n"
source "$(conda info --base)/etc/profile.d/conda.sh"

#CREATE ENVIRONMENT
echo "Creating Conda environment: $QFORTE_CONDA_ENV\n"
conda create -n "$QFORTE_CONDA_ENV" -y python=3.8

#INSTALL REQUIRED PACKAGES
echo "Installing Conda packages...\n"
conda install -n "$QFORTE_CONDA_ENV" -y -c conda-forge psi4 cmake openblas libopenblas pytest

#TODO: ADD AN OPTION TO INSTALL CUDA
#conda install -n $QFORTE_CONDA_ENV -c nvidia cuda

#ACTIVATE CONDA ENV
echo "Activating Conda environment...\n"
conda activate $QFORTE_CONDA_ENV

#SET THE CMAKE PREFIX TO THE CONDA PREFIX
echo "Setting CMAKE_PREFIX_PATH...\n"
sed -i "s|set(CMAKE_PREFIX_PATH \".*\")|set(CMAKE_PREFIX_PATH \"$CONDA_PREFIX\")|" CMakeLists.txt

#UPDATE THE MINIMUM REQUIRED VERSION OF CMAKE 
echo "Patching pybind11 and fmt CMakeLists.txt to use a supported version of CMake..."
#qforte
sed -i 's/^cmake_minimum_required(.*$/cmake_minimum_required(VERSION 3.5)/' CMakeLists.txt
#pybind11
sed -i 's/^cmake_minimum_required(.*$/cmake_minimum_required(VERSION 3.5)/' lib/pybind11/CMakeLists.txt
#fmt
sed -i 's/^cmake_minimum_required(.*$/cmake_minimum_required(VERSION 3.5)/' lib/fmt/CMakeLists.txt

#SET LIBOPENBLAS PATH DEPENDING ON OS
OS="$(uname)"
if [ "$OS" = "Linux" ]; then
	echo "You are on Linux: using libopenblas.so"
	sed -i 's|set(OPENBLAS_EXE ".*")|set(OPENBLAS_EXE ${CMAKE_PREFIX_PATH}/lib/libopenblas.so)|' CMakeLists.txt
elif [ "$OS" = "Darwin" ]; then
	echo "You are on MacOS: using libopenblas.dylib"
	sed -i 's|set(OPENBLAS_EXE ".*")|set(OPENBLAS_EXE ${CMAKE_PREFIX_PATH}/lib/libopenblas.dylib)|' CMakeLists.txt
	exit
else
	echo "Unknown OS -> exiting..."
	exit
fi

#BUILD
python setup.py develop && echo "qforte successfully installed" || { echo "qforte failed to install"; exit 1; }
