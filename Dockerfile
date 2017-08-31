#################################################################################
#
# DIRAC (testing) Docker for buildups with CI connected to private GitLab rep
#
# Written by Miro Ilias, GSI/Darmstadt, August-September,2017
#
#################################################################################
#
#
# we have LTS Ubuntu
FROM ubuntu:16.04
#
RUN apt-get update 
RUN apt-get install -y sudo 
RUN rm -rf /var/lib/apt/lists/* && sudo rm /var/cache/apt/archives/lock 
RUN sudo dpkg --configure -a
RUN sudo apt-get install cmake make
RUN sudo apt-get install gcc g++ gfortran
RUN sudo apt-get install zlib1g-dev 
RUN sudo apt-get remove openmpi-bin libopenmpi-dev 
RUN sudo apt-get install mpich2
# clean packages
RUN sudo apt-get autoclean
RUN sudo apt-get clean
RUN sudo apt-get autoremove
