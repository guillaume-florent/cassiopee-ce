#!/usr/bin/env bash

# File
#     install_cassiopee.sh
#
# Description
#     this script will
#     1) Build the cassiopee:2.5 Docker image
#     2) Create a container with the name cassiopee
#
#     Note: Docker daemon should be running before  launching script


username="$USER"
user="$(id -u)"
home="${1:-$HOME}"

imageName="guillaume-florent/cassiopee:2.5u1604"
containerName="cassiopee"
displayVar="$DISPLAY"

echo "*******************************************************"
echo ""
echo "Building Docker XFoil container ${containerName}"

docker build --tag ${imageName} .

docker run  -it -d --name ${containerName}                  \
    -e DISPLAY=${displayVar}                                \
    --workdir="${home}"                                     \
    --volume="${home}:${home}"                              \
     -v=/tmp/.X11-unix:/tmp/.X11-unix ${imageName}          \
     /bin/bash


echo "Container ${containerName} was created."

echo "********************************************************"
echo "Run the ./start_cassiopee.sh script to launch container "
echo "********************************************************"