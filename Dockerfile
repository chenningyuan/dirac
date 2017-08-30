FROM ubuntu:16.04
RUN apt-get update && apt-get install -y sudo && rm -rf /var/lib/apt/lists/* && sudo rm /var/cache/apt/archives/lock && sudo dpkg --configure -a
RUN sudo apt-get install cmake gcc g++ gfortran make zlib1g-dev & sudo apt-get remove openmpi-bin libopenmpi-dev & sudo apt-get install mpich2
