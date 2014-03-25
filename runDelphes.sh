#!/bin/bash
set -x

WORKDIR=${PWD}
LHDIR="/afs/cern.ch/work/o/obondu/LesHouches2013/"
DELPHES="Delphes-3.0.12"
CMSSW_release="CMSSW_5_3_9"
INPUT_FILE="/eos/cms/store/cmst3/user/obondu/HH/signal/13TeV/MGraviton_900_HHtobbbb.hepmc.gz"
#CONFIGURATION_FILE="CMS_cards/CMS_Phase_I_NoPileUp.tcl"
CONFIGURATION_FILE="CMS_cards/CMS_Phase_I_50PileUp.tcl"
OUTPUT_NAME="test"
ifile=0
OUTPUT_FILE="${OUTPUT_NAME}_${ifile}.root"


echo ${PWD}
ls

# Get what Delphes needs to run
cp -v  ${LHDIR}/${DELPHES}/DelphesHepMC .
cp -rv ${LHDIR}/${DELPHES}/CMS_cards/ .
cp -v  ${LHDIR}/${DELPHES}/MinBias.pileup .

# Setup ROOT and such via CMSSW
cd ${LHDIR}/${CMSSW_release}/src
eval `scram runtime -sh`
cd ${WORKDIR}

echo ${PWD}
ls

# Run the actual thing
xrd eoscms cat ${INPUT_FILE} | gunzip | ./DelphesHepMC ${CONFIGURATION_FILE} ${OUTPUT_FILE}

# Copy the stuff on eos
# TODO



