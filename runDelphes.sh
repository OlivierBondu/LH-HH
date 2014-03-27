#!/bin/bash

LHDIR="/afs/cern.ch/work/o/obondu/LesHouches2013/"
DELPHES="Delphes-3.0.12"
CMSSW_release="CMSSW_5_3_9"
INPUT_FILE="/eos/cms/store/cmst3/user/obondu/HH/signal/13TeV/MGraviton_900_HHtobbbb.hepmc.gz"
CONFIGURATION_FILE="Snowmass_cards/delphes_card_Snowmass_50PileUp.tcl"
PILEUP_FILE="Snowmass_pileup/MinBias100K_13TeV.pileup"
OUTPUT_NAME="test"
ifile=0

while getopts l:d:C:i:c:p:o:v: opt 
do
    case "${opt}" in 
    #l is LHDIR
        l)
            LHDIR=${OPTARG}
            ;;
    #d is DELPHES
        d)
            DELPHES=${OPTARG}
            ;;
    #C is CMSSW
        C)
            CMSSW_release=${OPTARG}
            ;;
    #i is input file
        i)
            INPUT_FILE=${OPTARG}
            ;;
    #c is config file
        c)
            CONFIGURATION_FILE=${OPTARG}
            ;;
    #p is pileup file
        p)
            PILEUP_FILE=${OPTARG}
            ;;
    #o is output name
        o)
            OUTPUT_NAME=${OPTARG}
            ;;
    #v is version
        v)
            ifile=${OPTARG}
            ;;
    # Error handling
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

WORKDIR=${PWD}
OUTPUT_FILE="${OUTPUT_NAME}_${ifile}.root"
echo "LHDIR= ${LHDIR}"
echo "DELPHES= ${DELPHES}"
echo "CMSSW_release= ${CMSSW_release}"
echo "INPUT_FILE= ${INPUT_FILE}"
echo "CONFIGURATION_FILE= ${CONFIGURATION_FILE}"
echo "PILEUP_FILE= ${PILEUP_FILE}"
echo "OUTPUT_NAME= ${OUTPUT_NAME}"


echo ${PWD}
ls

# Get what Delphes needs to run
cp -v  ${LHDIR}/${DELPHES}/DelphesHepMC .
cp -rv ${LHDIR}/${DELPHES}/Snowmass_cards/ .
cp -v  ${LHDIR}/${DELPHES}/${PILEUP_FILE} ./MinBias.pileup

# Setup ROOT and such via CMSSW
cd ${LHDIR}/${CMSSW_release}/src
eval `scram runtime -sh`
cd ${WORKDIR}

echo ${PWD}
ls

set -x
# Run the actual thing
xrd eoscms cat ${INPUT_FILE} | gunzip | ./DelphesHepMC ${CONFIGURATION_FILE} ${OUTPUT_FILE}

# Copy the stuff on eos
# TODO
# cp -v ${OUTPUT_FILE} ${EOS_OUTPUT}



