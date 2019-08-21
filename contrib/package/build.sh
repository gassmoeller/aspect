#!/bin/bash

apt-get update && apt-get upgrade -yq && apt-get install -yq --no-install-recommends build-essential cmake gcc g++ gfortran git libblas-dev liblapack-dev libopenmpi-dev openmpi-bin numdiff zlib1g-dev fuse wget ca-certificates lsb-release file # libmpich-dev mpich

git clone https://github.com/dealii/candi.git/
cd candi
./candi.sh -j6 -p ${PWD}/../AppImage --packages="hdf5 p4est trilinos dealii"
cd ../AppImage
git clone https://github.com/geodynamics/aspect.git/
mkdir build && cd build
cmake -D ASPECT_PRECOMPILE_HEADERS=ON -DDEAL_II_DIR=../deal.II-v9.1.1 ../aspect
make release
make -j6

cd ../
wget https://github.com/linuxdeploy/linuxdeploy/releases/download/continuous/linuxdeploy-x86_64.AppImage
chmod +x linuxdeploy-x86_64.AppImage
./linuxdeploy-x86_64.AppImage --appimage-extract

export VERSION=dev
squashfs-root/AppRun -e ${PWD}/build/aspect -d ${PWD}/aspect/contrib/package/aspect.desktop -i ${PWD}/aspect/contrib/package/mesh-2d-icon.png --appdir AppDir --output appimage
squashfs-root/AppRun --appdir AppDir --output appimage
cp ASPECT-x86_64.AppImage software/
