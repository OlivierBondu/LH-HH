########################################
# SET ENVIRONMENT VARs FOR ROOT
export ROOTSYS=/home/amassiro/root-v5-34/
export PATH=$PATH:$ROOTSYS/bin
export PATH=$PATH:$ROOTSYS/include
# export LD_LIBRARY_PATH=$ROOTSYS/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib
########################################

##################################################
# new SET ENVIRONMENT VARs FOR ROOT ---> ROOT 5.25
. $ROOTSYS/bin/thisroot.sh
##################################################



export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:.