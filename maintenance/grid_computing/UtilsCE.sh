############################################################################
#
# List of specific function in shell language related to execution of DIRAC 
# software in grid environment.
#
# Add it to your script at the beginning as "source UtilsCE.sh"
# as function declarations must precede their calls.
#
# Written by Miroslav Ilias, Oct.2012
#
############################################################################

get_nprocs_CE()
{ # function , http://www.linuxjournal.com/content/return-values-bash-functions
  local result=$1 # store variable name
  unser nprocs
  local nprocs=`cat /proc/cpuinfo | grep processor | wc -l`
  eval $result="'$nprocs'"
}

print_CE_info ()
{ # prints basic info about given CE
unset THIS_HOSTNAME
THIS_HOSTNAME=`hostname -f`
if [[ -z $THIS_HOSTNAME ]]; then
 THIS_HOSTNAME=`hostname`
fi
echo -e "\n --- Hostname: $THIS_HOSTNAME "; 
echo -e "Server specification: \c";  uname -a

######################################################
#     get number of CPUs on the given node
######################################################
get_nprocs_CE nprocs_all

if [[ -n $PBS_NODEFILE ]]; then
nprocs_PBS=`cat $PBS_NODEFILE | wc -l`
fi
nprocs=$nprocs_all
echo -e "\n==== CPU \c"; cat /proc/cpuinfo | grep "model name"  | tail -n 1
echo -e "==== Number of all CPUs     :  $nprocs   ";
echo -e "==== Number of allowed CPUs :  $nprocs_PBS  ";
if [[ -n $PBS_NODEFILE ]]; then
  echo -e "PBS_NODEFILE=$PBS_NODEFILE"
  echo -e "PBS_NUM_PPN=$PBS_NUM_PPN"
fi
#####################
###  memory of CE ###
#####################
echo -e "\n ---- server's total memory :\c"; cat /proc/meminfo | grep MemTotal
echo -e "through the 'free -t -m' command (values in MB): "; free -m -t
#####################
# check space on CE #
#####################
echo -e "\n--- Server's filesystem (df): "; df -h .
echo -e "--- My current directory :  $PWD "
echo -e "--- Occupied space by me (du) : \c"; du -h -s
echo -e "\n\n ****  server environment:  set command  **** \n" 
set
echo -e "\n\n Is there cmake ? which cmake="; which cmake
echo -e "Python present : \c"; which python; python -V

} # Function declaration must precede call.

querry()
{
echo "..."
if [[ "$1" == "voce" ]]; then
  echo -e "\n Checking long-term proxy with myproxy-info -v : "
  myproxy-info -v -d --pshost px.ui.savba.sk
fi
}

querry_CE_attributes()
{
echo "..."
if [[ "$1" == "voce" ]]; then
  echo -e "\n Checking long-term proxy with myproxy-info -v : "
  myproxy-info -v -d --pshost px.ui.savba.sk
fi
if [[ "$1" == "compchem" ]]; then
  echo -e "\n Checking long-term proxy with myproxy-info -v : "
  myproxy-info -v -d --pshost myproxy.cnaf.infn.it
fi
echo -e "\n Checking short-term proxy with voms-proxy-info --all : "
voms-proxy-info --all
if [[ -n $CE_ID ]]; then
# get info about CE #
echo -e "\n querry attribute SMPSize, CE=$CE_ID:"
date
lcg-info --list-ce --query "CE=$CE_ID" --attrs SMPSize
date
echo -e "\n querry attribute PlatformArch:"
lcg-info --list-ce --query "CE=$CE_ID" --attrs PlatformArch
echo -e "querry attributes TotalCPUs, Memory, ClockSpeed, MaxCPUTime :"
lcg-info --list-ce --query "CE=$CE_ID" --attrs TotalCPUs,Memory,ClockSpeed,MaxCPUTime
echo -e "\n faster querrying using CE_BDII=$CE_BDII....";date
CE_BDII=${CE_ID/:*/:2170}
lcg-info --list-ce --bdii "$CE_BDII" --query "CE=$CE_ID" --attrs SMPSize,TotalCPUs,PlatformArch,Memory,ClockSpeed,MaxCPUTime
date
fi
}

check_file_on_SE()
{ # parametric call with name of VO
VO=$1
echo -e "\n\n On the host CE, checking variable LFC_HOST=$LFC_HOST"
if [[ -z $LFC_HOST ]];  then
  echo -e "The LFC_HOST environmental variable is not defined !"
fi
echo -e "Trying command lcg-infosites --vo $VO lfc :"; lcg-infosites --vo $VO lfc

echo -e  "\n Files in my vo=$VO SE space, lcg-ls -l lfn://grid/$VO/ilias:"
lcg-ls -l lfn://grid/$VO/ilias || { echo -e "...previous command failed, once again... "; lcg-ls -l lfn://grid/$VO/ilias; }
echo -e "\n  File attributes lcg-lr  lfn://grid/$VO/ilias/dirac_grid_suite.tgz:"
lcg-lr    lfn://grid/$VO/ilias/dirac_grid_suite.tgz || { echo -e "...previous command failed, once again... "; lcg-lr  lfn://grid/$VO/ilias/dirac_grid_suite.tgz;  }
echo " ------------------------------------------------------ "
lfc-getacl  /grid/$VO/ilias
}

download_from_SE()
{ # download file to the current CE, parameter - name of VO
VO=$1
echo -e "\n Going to download from SE, lcg-cp lfn://grid/$VO/ilias/dirac_grid_suite.tgz  file://$PWD/dirac_grid_suite.tgz...."
lcg-cp -v lfn://grid/$VO/ilias/dirac_grid_suite.tgz  file://$PWD/dirac_grid_suite.tgz || { echo "Error exit,downloading of dirac_grid_suite.tgz failed - once again" ; lcg-cp -v lfn://grid/$VO/ilias/dirac_grid_suite.tgz  file://$PWD/dirac_grid_suite.tgz || {echo -e "Download failed second time ! Exit !!! "; exit -4} }

echo -e "\nDownloaded file attributes (ls -lt dirac_grid_suite.tgz) :"
ls -lt dirac_grid_suite.tgz || { echo -e "\n\nNo dowloaded dirac_grid_suite.tgz file found in CE current directory - error exit !"; exit -3;}
echo " ----------------------------------------------------------------------------------------------------"
}

final_message()
{ # flush out final message
echo -e "\n -------------------------------------------------------"
echo -e "Job finished at: \c"; date
echo -e "I am here: \c"; pwd
echo -e "My space: "; df -h .
echo -e "\n All files:"; ls -lt
echo -e "\n Total occupied space: \c"; du -h -s .
}
