#!/bin/bash

# Use cmake3.28  or greater
#export PATH=$HOME/packages/cmake-3.28.3/bin:$PATH

export BUILD_HOME=`pwd`
echo "BUILD_HOME: $BUILD_HOME"

. /net/projects/tools/x86_64/rhel-8/intel-oneapi/2023.1/setvars.sh > /dev/null

export PATH=/usr/local/cuda-11.7/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda-11.7/lib64:$LD_LIBRARY_PATH

cd $BUILD_HOME
rm -rf ./umt_workspace/build_conduit
rm ./umt_workspace/CMakeCache.txt
rm -rf ./umt_workspace/teton/CMakeFiles
cd umt_workspace
pwd
#make clean
cd $BUILD_HOME

echo "Building ..."
./build_and_run_umt.sh > ./build.log 2>&1
