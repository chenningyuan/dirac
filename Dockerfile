FROM ubuntu:16.04
RUN apt-get update && apt-get install -y sudo && rm -rf /var/lib/apt/lists/*
RUN sudo apt-get install gcc g++ gfortran make zlib1g-dev & sudo apt-get remove openmpi-bin libopenmpi-dev & sudo apt-get install mpich2
