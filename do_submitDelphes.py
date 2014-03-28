#!/usr/bin/env python
import subprocess

queue = '8nh'
eos = "/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select"
eos_base = "/eos/cms/store/cmst3/user/obondu/HH/delphes/"
# version
version = 0
# center of mass energies
com = "13TeV"
# PU scenario
#PU = 'NoPileUp'
#PU = '50PileUp'
PU = '140PileUp'
# samples (final state, mass hyp)
eos_input = "/store/cmst3/user/obondu/HH/signal/13TeV/"
list_of_files = subprocess.check_output([eos, "ls", eos_input])
list_of_files = [x for x in list_of_files.split('\n') if "hepmc" in x and "bbbb" in x]
#sample_files = [ eos_input + x for x in list_of_files ]
#sample_files = [ "/store/cmst3/user/obondu/HH/signal/13TeV/MGraviton_900_HHtoWWbb.hepmc.gz" ]

tosubmit = []
# Create the list of submission commands
for file in list_of_files:
    sample = eos_input + file
    outname = file.replace(".hepmc.gz", "")
    eos_output = eos_base + "/" + str(com) + "/" + str(PU) + "/"
    tosubmit.append("bsub -q " + str(queue) + "\trunDelphes.sh" + "\t-i " + str(sample) + "\t-c Snowmass_cards/delphes_card_Snowmass_" + str(PU) + ".tcl" + "\t-p Snowmass_pileup/MinBias100K_" + str(com) + ".pileup" + "\t-O " + str(eos_output) + "\t-o " + str(outname) + "\t-v " + str(version))

# Execute the submission list
for job in tosubmit:
    print job
    subprocess.Popen(job, stdout=subprocess.PIPE, shell=True)
