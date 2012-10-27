############################################################################
#
# List of specific functions in shell language related to execution of DIRAC 
# software in grid environment.
#
# Add it to your script at the beginning as "source UtilsCE.sh"
# as function declarations must precede their calls.
#
# Written by Miroslav Ilias, Oct.2012
#
############################################################################

echo -e "\n Going to source - include - file with bash shell functions, $0 ."

function get_nprocs_CE()
{ # function , http://www.linuxjournal.com/content/return-values-bash-functions
  local result=$1 # store variable name
  unset nprocs
  local nprocs=`cat /proc/cpuinfo | grep processor | wc -l`
  eval $result="'$nprocs'"
}

function print_CE_info ()
{ # prints basic info about given CE
unset THIS_HOSTNAME
THIS_HOSTNAME=`hostname -f`
if [[ -z $THIS_HOSTNAME ]]; then
 THIS_HOSTNAME=`hostname`
fi
echo -e "\n --- Hostname: $THIS_HOSTNAME "; 
echo -e "Server specification: \c";  uname -a
echo -e "\n shell version, sh --version:"; sh --version
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
return 0
} 

function querry_CE_attributes()
{ #   ask on values of CE parameters   # 
local VO=$1
echo -e "\n First, on CE=$CE, VO=$VO checking long-term proxy with myproxy-info -v : "
if [[ "$VO" == "voce" ]]; then
  myproxy-info -v -d --pshost px.ui.savba.sk
elif [[ "$VO" == "compchem" ]]; then
  myproxy-info -v -d --pshost myproxy.cnaf.infn.it
else
  echo -e "\n What is your VO ? Neither voce,nor compchem..."
fi
#
echo -e "\n Checking short-term proxy, voms-proxy-info --all : "
voms-proxy-info --all
#
if [[ -n $CE_ID ]]; then
# get info about CE #
  echo -e "\n querrying attribute on this CE=$CE_ID:"
  date
  echo -e "querry attributes, like  PlatformArch, TotalCPUs, Memory, ClockSpeed, MaxCPUTime...."
  lcg-info --list-ce --query "CE=$CE_ID" --attrs PlatformArch,SMPSize,TotalCPUs,Memory,VMemory,ClockSpeed,MaxCPUTime
  CE_BDII=${CE_ID/:*/:2170}
  echo -e "\n Now do faster querrying using CE_BDII=$CE_BDII";date
  lcg-info --list-ce --bdii "$CE_BDII" --query "CE=$CE_ID" --attrs SMPSize,TotalCPUs,PlatformArch,Memory,VMemory,ClockSpeed,MaxCPUTime,CPUVendor,PhysicalCPU,CE,CEStatus,FreeCPUs,Cluster,VOInfoVOs
  date
else
  echo -e "\n No querrying of attributes as CE_ID not defined !"
fi
return 0
}

function download_from_SE()
{ # download file from SE to the current CE local directory, parameter - name of VO
local VO=$1
#local getcmd="lcg-cp -v lfn://grid/$VO/ilias/dirac_grid_suite.tgz file://$PWD/dirac_grid_suite.tgz"
# remove verbose mode to eliminate stderr streams which are occuring 
local getcmd="lcg-cp lfn://grid/$VO/ilias/dirac_grid_suite.tgz file://$PWD/dirac_grid_suite.tgz"
echo -e "\n Going to download from SE, 1st attempt, $getcmd "
$getcmd
RETVAL=$?
# this does not work..each time it's repeating download..
#RETVAL=$?; [ $RETVAL -ne 0 ] && echo "some error, downloading once again ";$getcmd 
#RETVAL=$?; [ $RETVAL -ne 0 ] && echo "downloading failed - exit now "; exit 17
echo -e "FYI: downloading ended with code $RETVAL"

# this command does not work in bash - too compicated !
#$getcmd || { echo "once again"; {$getcm || {echo -e "Download failed second time ! Exit !!! "; exit 41} }}
#
echo -e "\nDownloaded file attributes (ls -lt dirac_grid_suite.tgz) :"
ls -lt dirac_grid_suite.tgz || { echo -e "\n\nNo dowloaded file found in CE current directory - exit !"; exit 32;}
echo " -----------------------------------------------------------------------------------------"
return 0
}

function check_file_on_SE()
{ # parametric call with name of VO
local VO=$1
echo -e "\n\n On the host CE, checking variable LFC_HOST=$LFC_HOST"
if [[ -z $LFC_HOST ]];  then
  echo -e "The LFC_HOST environmental variable is not defined !"
  echo -e "Trying command lcg-infosites --vo $VO lfc :"; lcg-infosites --vo $VO lfc
  export LFC_HOST=`lcg-infosites --vo $VO lfc`
fi
echo -e  "\n Files in my vo=$VO SE space, lcg-ls -l lfn://grid/$VO/ilias:"
lcg-ls -l lfn://grid/$VO/ilias || { echo -e "...previous command failed, once again... "; lcg-ls -l lfn://grid/$VO/ilias; }
echo -e "\n File attributes lcg-lr lfn://grid/$VO/ilias/dirac_grid_suite.tgz:"
lcg-lr lfn://grid/$VO/ilias/dirac_grid_suite.tgz || { echo -e "...previous command failed, once again... "; lcg-lr  lfn://grid/$VO/ilias/dirac_grid_suite.tgz;  }
echo " ------------------------------------------------------ "
# for these commands we need the LFC_HOST variable
echo -e "\n check users's LFC-directory, lfc-getacl  /grid/$VO/ilias:"
lfc-getacl  /grid/$VO/ilias
echo -e "\n check user's LFC-file,  lfc-getacl  /grid/$VO/ilias/dirac_grid_suite.tgz:"
lfc-getacl  /grid/$VO/ilias/dirac_grid_suite.tgz
echo -e "-----------------------------------------------------------\n"
return 0
}

function unpack_DIRAC()
{ # unpack and check the freshly downloaded DIRAC tarball
echo -e "\n\n ---- Unpacking dowloaded DIRAC suite, tar : \c "; which tar;  tar --version
echo -e "performing tar xfz dirac_grid_suite.tgz ...\c"
tar xfz dirac_grid_suite.tgz || {echo "problem in unpacking!"; exit 12}
echo -e "unpacking done; see the content of $PWD :"; ls -lt
#echo -e "\n making executables:"
#ls -lt $BUILD/dirac.x  $BUILD/cmake  $BUILD/ctest $BUILD/mpirun  pam  runtest pamadm setup
#chmod  u+x  $BUILD/dirac.x  $BUILD/cmake  $BUILD/ctest $BUILD/mpirun  pam  runtest pamadm setup
#ls -lt $BUILD/dirac.x  $BUILD/cmake  $BUILD/ctest $BUILD/mpirun  pam  runtest pamadm setup
# create diracrc file for parallel runs on the grid CE
cp maintenance/grid_computing/diracrc_grid  diracrc
echo -e "\n ---- diracrc created: "; cat diracrc
return 0
}

function final_message()
{ 
echo -e "\n -------------------------------------------------------"
echo -e "Job finished at: \c"; date
echo -e "Filesystem on my space: "; df -h $PWD
echo -e "\n All files:"; ls -lt
echo -e "\n Total occupied space: \c"; du -h -s $PWD
echo -e "\n ------------------------------------------------------------------"
return 0
}

echo -e "List of declated bash shell functions in sourced file $0 :"
declare -F
