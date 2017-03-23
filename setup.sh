#!/bin/bash

if [ $KP_ANALYZER_PATH ]; then
    echo "KP_ANALYZER_PATH is already defined: use a clean shell!"
    return 1
fi

export KP_ANALYZER_PATH=$(pwd)

if [ $HOSTNAME == "charm" ]; then
	source /usr/local/Geant4/bin/geant4.sh
	#scl enable devtoolset-2 bash # -- enable gcc 4.8 -- #
fi

echo "================ environment ================"
echo "KP_ANALYZER_PATH:" $KP_ANALYZER_PATH
echo "G4NEUTRONHPDATA:" $G4NEUTRONHPDATA
echo "G4LEDATA:" $G4LEDATA
echo "G4LEVELGAMMADATA:" $G4LEVELGAMMADATA
echo "G4RADIOACTIVEDATA:" $G4RADIOACTIVEDATA
echo "G4PIIDATA:" $G4PIIDATA
echo "G4REALSURFACEDATA:" $G4REALSURFACEDATA
echo "G4SAIDXSDATA:" $G4SAIDXSDATA
echo "G4ABLADATA:" $G4ABLADATA
echo "G4ENSDFSTATEDATA:" $G4ENSDFSTATEDATA
if [ -z $G4ENSDFSTATEDATA ]; then
    echo "     [WARNING]: Geant4 is not available in this machine"
fi
echo "============================================="
echo "setup is finished. Welcome :)"

if [ $HOSTNAME == "charm" ]; then
	echo "Use gcc 4.8 at CentOS6: scl enable devtoolset-2 bash"
	scl enable devtoolset-2 bash # -- enable gcc 4.8 -- #
fi