FROM docker.io/debian:bookworm-20231030 as debian-base
MAINTAINER Martin Heistermann <martin.heistermann@inf.unibe.ch>

FROM debian-base as base-packages

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    binutils \
    build-essential \
    ca-certificates \
    clang \
    cmake \
    curl \
    git \
    libc-dev \
    liblapack-dev \
    libopenblas64-serial-dev \
    libtool \
    locales \
    ninja-build \
    time \
    tzdata \
    wget \
    pkg-config \
    libgmp-dev \
    libomp-dev \
    gfortran \
    clang

RUN echo "Downloading Gurobi..." && \
            mkdir -p /opt && \
            cd /opt && \
            wget -q https://packages.gurobi.com/10.0/gurobi10.0.3_linux64.tar.gz && \
            echo 82f916db110c42ce8ce13c10a14eba97c7acd63c3c0c59f98186c5085780ca83  gurobi10.0.3_linux64.tar.gz | sha256sum --check && \ 
            tar xf gurobi10.0.3_linux64.tar.gz && \
            ln -s gurobi1003 gurobi

RUN echo "Downloading Linuxdeploy..." && \
    wget -q https://github.com/linuxdeploy/linuxdeploy/releases/download/continuous/linuxdeploy-x86_64.AppImage -o /opt/linuxdeploy && chmod +x /opt/linuxdeploy

RUN echo "Downloading coinbrew (for ipopt)" && \
    mkdir -p /usr/src/coin-or && \
    mkdir -p /opt/coin-or && \
    cd /usr/src/coin-or && \
    wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew \
    && chmod +x coinbrew

RUN cd /usr/src/coin-or && ./coinbrew fetch Ipopt@3.14.13
RUN cd /usr/src/coin-or && ./coinbrew build Ipopt@3.14.13 --prefix=/opt/coin-or --parallel-jobs 8
# TODO: can set --parallel-jobs in build step

FROM base-packages AS build

COPY . /app

ENV IPOPT_HOME=/opt/coin-or

RUN mkdir /app/build && cd /app/build && \
    cmake -G Ninja \
    -D CMAKE_BUILD_TYPE=Release \
    -D CMAKE_CXX_COMPILER=clang++ \
    -D CMAKE_C_COMPILER=clang \
    -D CMAKE_CXX_FLAGS="-march=native" \
    ..
#-D CMAKE_INTERPROCEDURAL_OPTIMIZATION=ON \

RUN ninja
RUN ninja install DESTDIR=AppDir

RUN /opt/linuxdeploy --appdir AppDir --output appimage

RUN ln -s /app/build/Build/bin/* /usr/local/bin
#RUN apt-get clean && rm -rf /var/lib/apt/lists/

