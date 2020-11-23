FROM ubuntu:20.04
# Fix the installation of tzdata for Ubuntu 20.04
ARG TIMEZONE=Europe/Berlin
RUN export TZ=$TIMEZONE && echo $TZ > /etc/timezone && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    apt-get -yy update && apt-get -yy install wget tzdata

RUN apt-get -yy install build-essential cmake libeigen3-dev libxml2-dev libboost-all-dev petsc-dev python3-dev python3-numpy python3-pip git-all

RUN wget -q -O libprecice.deb https://github.com/precice/precice/releases/download/v2.1.1/libprecice2_2.1.1_focal.deb && \
    ( dpkg -i libprecice.deb || apt-get install -fyy ) && \
    rm libprecice.deb && \
    binprecice xml > /dev/null # Make sure the installation is functional
    
RUN pip3 install pyprecice --user

RUN pip3 install matplotlib nutils --user

RUN apt-get -yy install software-properties-common && add-apt-repository -yy ppa:fenics-packages/fenics && apt-get -yy update && apt-get -yy install --no-install-recommends fenics

RUN pip3 install fenicsprecice --user
