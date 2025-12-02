FROM docker.io/debian:trixie-20251117 as debian-base
MAINTAINER Martin Heistermann <martin.heistermann@unibe.ch>

FROM debian-base as base-packages

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    binutils \
    build-essential \
    ca-certificates \
    g++ \
    cmake \
    curl \
    git \
    libc-dev \
    liblapack-dev \
    libopenblas64-serial-dev \
    libopenmpi-dev \
    libtool \
    locales \
    ninja-build \
    time \
    tzdata \
    wget \
    pkg-config \
    libgmp-dev \
    libsuitesparse-dev \
    gfortran

# Currently (2025-12-02) unsuitable without hacks:
# coinor-libipopt-dev
# coinor-libbonmin-dev \

RUN mkdir -p /usr/src/coin-or && \
    mkdir -p /opt/coin-or && \
    cd /usr/src/coin-or && \
    wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew \
    && chmod +x coinbrew

RUN cd /usr/src/coin-or && ./coinbrew fetch https://github.com/coin-or-tools/ThirdParty-Mumps@3.0.11

RUN cd /usr/src/coin-or && ./coinbrew fetch Ipopt@3.14.19 --skip-update
RUN cd /usr/src/coin-or && ./coinbrew fetch Bonmin@master --skip-update

RUN cd /usr/src/coin-or && ./coinbrew build Ipopt@3.14.19 --verbosity 2 --skip-update --prefix=/opt/coin-or --parallel-jobs 8 --tests none
RUN cd /usr/src/coin-or && ./coinbrew build Bonmin@master --verbosity 2 --skip-update --prefix=/opt/coin-or --parallel-jobs 8 --tests none

RUN ln -s /opt/coin-or/include/coin-or /opt/coin-or/include/coin

FROM base-packages AS build

COPY . /app

ENV CoinUtils_DIR=/opt/coin-or
ENV IPOPT_HOME=/opt/coin-or
ENV CBC_DIR=/opt/coin-or
ENV CLP_DIR=/opt/coin-or

RUN mkdir /app/build && cd /app/build && \
    cmake -G Ninja \
    -D CMAKE_BUILD_TYPE=RelWithDebInfo \
    -D BONMIN_ROOT_DIR=/opt/coin-or \
    ..
#-D CMAKE_INTERPROCEDURAL_OPTIMIZATION=ON \

RUN cd /app/build && ninja

RUN ln -s /app/build/Build/bin/* /usr/local/bin
#RUN apt-get clean && rm -rf /var/lib/apt/lists/

