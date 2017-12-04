#!/usr/bin/env bash

# File
#     start_cassiopee.sh
#
# Description
#      This script will
#          1) Start the cassiopee container with name 'cassiopee'
#             in the the shell terminal.
#      Note
#          1) User should run xhost + from other terminal
#          2) Docker daemon should be running before launching this script
#
# At the container prompt :
# > source env_Cassiopee.sh (in /opt/cassiopee/build/Dist)
# > source /opt/cassiopee/build/Dist/env_Cassiopee.sh
# > python
#   Test installation by trying to import KCore
# OR
# run cassiopee (in /opt/cassiopee/build/Dist/bin/x86_r8)
#
#------------------------------------------------------------------------------

xhost +local:cassiopee
docker start cassiopee
docker exec -it cassiopee /bin/bash
# or
# docker exec -it cassiopee cassiopee