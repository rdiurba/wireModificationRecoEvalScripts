#!/bin/bash

STEM=$1
macrofile=$2 
inputfile=$3 

# Create a dummy file in the output directory. This is a hack to get jobs
# that fail to end themselves quickly rather than hanging on for a long time
# waiting for output to arrive.
DUMMY_OUTPUT_FILE=${CONDOR_DIR_OUTPUT}/${JOBSUBJOBID}_${STEM}_dummy_output
touch ${DUMMY_OUTPUT_FILE}

# Get the source code for GENIE v3
cd $CONDOR_DIR_INPUT
source /cvmfs/larsoft.opensciencegrid.org/products/setups
setup git v2_15_1

setup root v6_12_04e -q e15:prof
setup lhapdf v5_9_1k -q e15:prof
setup pdfsets v5_9_1b
setup log4cpp v1_1_3a -q e15:prof


export _RUN_NUMBER=$((PROCESS+1))
echo "Running $0 on "$HOSTNAME
echo "Cluster: " ${CLUSTER}
echo "Process: " ${PROCESS}

echo " Sourcing everything...."

source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh

setup uboonecode v08_00_00_19 -q e17:prof
setup fife_utils
echo "Dumping environment."

touch ${DUMMY_OUTPUT_FILE}
#mv ${_RUN_NUMBER}.out ${_CONDOR_SCRATCH_DIR}/.
#cd ${_CONDOR_SCRATCH_DIR}

echo ${_RUN_NUMBER}
echo ${macroFile}
echo ${inputfile}
echo "fullpath=`head -n ${_RUN_NUMBER} ${inputfile} | tail -n 1`"
fullpath=`head -n ${_RUN_NUMBER} ${inputfile}  | tail -n 1`
echo $fullpath

echo " Copying File...."
ifdh cp $fullpath ${PWD}/test.root
ls ${PWD}

echo "What's in here? "
ls

echo " Run time! "

root -l -b -q ${macrofile} 

mkdir ${_RUN_NUMBER}

#cp ${_RUN_NUMBER}.out ${_RUN_NUMBER}/.
cp output.root ${_RUN_NUMBER}/${_RUN_NUMBER}.root

ifdh mkdir ${CONDOR_DIR_OUTPUT}
ifdh mkdir ${CONDOR_DIR_OUTPUT}/${_RUN_NUMBER}
ifdh cp -r ${_RUN_NUMBER} ${CONDOR_DIR_OUTPUT}/



