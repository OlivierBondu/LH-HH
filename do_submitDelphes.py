#!/usr/bin/env python
import subprocess

queue = '8nh'

eos_base = "/eos/cms/store/cmst3/user/obondu/HH/delphes/"
# version
version = 0
# center of mass energies
com = "13TeV"
# PU scenario
PU = 'NoPileUp'
# samples (final state, mass hyp)
sample_files = [ "/store/cmst3/user/obondu/HH/signal/13TeV/MGraviton_900_HHtoWWbb.hepmc.gz" ]

tosubmit = []
# Create the list of submission commands
for sample in sample_files:
    tosubmit.append("bsub -q " + str(queue) + "\trunDelphes.sh" + "\t-i " + str(sample) + "\t-c Snowmass_cards/delphes_card_Snowmass_" + str(PU) + ".tcl" + "\t-p Snowmass_pileup/MinBias100K_" + str(com) + ".pileup" + "\t-O " + str(eos_base) + "\t-o " + "TEST" + "\t-v " + str(version))

# Execute the submission list
for job in tosubmit:
    print job
    subprocess.Popen(job, stdout=subprocess.PIPE, shell=True)

