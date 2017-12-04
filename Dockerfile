FROM ubuntu:16.04

MAINTAINER Guillaume Florent version: 0.1

RUN apt-get update && apt-get install -y \
    g++ \
    git \
    libgl1-mesa-dev \
    mesa-common-dev \
    mesa-utils \
    freeglut3 \
    freeglut3-dev \
    libglew1.5 \
    libglew1.5-dev \
    libglu1-mesa \
    libglu1-mesa-dev \
    libgl1-mesa-glx \
    python-dev \
    python-tk \
    python-numpy \
    scons \
    gfortran \
    xorg-dev \
    && rm -rf /var/lib/apt/lists/*

ENV CASSIOPEE /opt/cassiopee-ce/build
ENV LD_LIBRARY_PATH /usr/lib/x86_64-linux-gnu

WORKDIR /opt
ADD https://api.github.com/repos/guillaume-florent/cassiopee-ce/git/refs/heads/master version.json
RUN git clone https://github.com/guillaume-florent/cassiopee-ce

WORKDIR /opt/cassiopee-ce/cassiopee/KCore/
RUN ./install && . $CASSIOPEE/Dist/env_Cassiopee.sh && ./install

WORKDIR /opt/cassiopee-ce/cassiopee/Compressor/
RUN ./install && . $CASSIOPEE/Dist/env_Cassiopee.sh && ./install

WORKDIR /opt/cassiopee-ce/cassiopee/Connector/
RUN ./install && . $CASSIOPEE/Dist/env_Cassiopee.sh && ./install

WORKDIR /opt/cassiopee-ce/cassiopee/Converter/
RUN ./install && . $CASSIOPEE/Dist/env_Cassiopee.sh && ./install

WORKDIR /opt/cassiopee-ce/cassiopee/Dist2Walls/
RUN ./install && . $CASSIOPEE/Dist/env_Cassiopee.sh && ./install

WORKDIR /opt/cassiopee-ce/cassiopee/Distributor2/
RUN ./install && . $CASSIOPEE/Dist/env_Cassiopee.sh && ./install

WORKDIR /opt/cassiopee-ce/cassiopee/Generator/
RUN ./install && . $CASSIOPEE/Dist/env_Cassiopee.sh && ./install

WORKDIR /opt/cassiopee-ce/cassiopee/Geom/
RUN ./install && . $CASSIOPEE/Dist/env_Cassiopee.sh && ./install

WORKDIR /opt/cassiopee-ce/cassiopee/Initiator/
RUN ./install && . $CASSIOPEE/Dist/env_Cassiopee.sh && ./install

WORKDIR /opt/cassiopee-ce/cassiopee/Intersector/
RUN ./install && . $CASSIOPEE/Dist/env_Cassiopee.sh && ./install

WORKDIR /opt/cassiopee-ce/cassiopee/Post/
RUN ./install && . $CASSIOPEE/Dist/env_Cassiopee.sh && ./install

WORKDIR /opt/cassiopee-ce/cassiopee/RigidMotion/
RUN ./install && . $CASSIOPEE/Dist/env_Cassiopee.sh && ./install

WORKDIR /opt/cassiopee-ce/cassiopee/Transform/
RUN ./install && . $CASSIOPEE/Dist/env_Cassiopee.sh && ./install

WORKDIR /opt/cassiopee-ce/cassiopee/CPlot/
RUN ./install && . $CASSIOPEE/Dist/env_Cassiopee.sh && ./install

WORKDIR /opt/cassiopee-ce/build/Dist/

CMD ["/bin/bash"]
