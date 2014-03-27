LH-HH
=====
  
  Two things needs to be installed: Delphes and this package. Delphes needs CMSSW to be run (see below), we'll assume every of this three things is located in the same directory.

## First: run Delphes

### Setting things up

Delphes needs ROOT to run, we'll just go through CMSSW for that

        cmsrel CMSSW_5_3_9
        cd CMSSW_5_3_9/src
        cmsenv
        cd ../..
        
Next you need to get / compile Delphes 3.0.12 (as of 2013-03-25)

        wget http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.0.12.tar.gz
        tar -xvzf Delphes-3.0.12.tar.gz
        cd Delphes-3.0.12
        make -j 8

Once Delphes is setup, we want to get some realistic detector parameters. Some sets can be found from the Snowmass exercise 2013 [here](http://www.snowmass2013.org/tiki-index.php?page=Energy_Frontier_FastSimulation) and referencable from [arXiv:1309.1057](http://arxiv.org/abs/1309.1057). Below how to get the cards themselves

        # Get detector configuration cards
        mkdir Snowmass_cards
        wget -O Snowmass_cards/delphes_card_Snowmass_NoPileUp.tcl http://cvs.web.cern.ch/cvs/cgi-bin/viewcvs.cgi/UserCode/spadhi/Snowmass/Cards/delphes_card_Snowmass_NoPileUp.tcl?view=co
        wget -O Snowmass_cards/delphes_card_Snowmass_50PileUp.tcl http://cvs.web.cern.ch/cvs/cgi-bin/viewcvs.cgi/UserCode/spadhi/Snowmass/Cards/delphes_card_Snowmass_50PileUp.tcl?view=co
        wget -O Snowmass_cards/delphes_card_Snowmass_140PileUp.tcl http://cvs.web.cern.ch/cvs/cgi-bin/viewcvs.cgi/UserCode/spadhi/Snowmass/Cards/delphes_card_Snowmass_140PileUp.tcl?view=co
        wget -O Snowmass_cards/delphes_card_Snowmass_VLHCPileUp.tcl http://cvs.web.cern.ch/cvs/cgi-bin/viewcvs.cgi/UserCode/spadhi/Snowmass/Cards/delphes_card_Snowmass_VLHCPileUp.tcl?view=co
        # Correct some 3.0.10 to 3.0.12 Delphes versioning
        sed -i -e 's/AdjacencyCut 2\.0/AdjacencyCut 2/g' Snowmass_cards/delphes_card_Snowmass_NoPileUp.tcl
        # Get pileup files
        mkdir Snowmass_pileup
        wget -P Snowmass_pileup http://red-gridftp11.unl.edu/Snowmass/MinBias100K_14TeV.pileup
        wget -P Snowmass_pileup http://red-gridftp11.unl.edu/Snowmass/MinBias100K_13TeV.pileup
        wget -P Snowmass_pileup http://red-gridftp11.unl.edu/Snowmass/MinBias100K_33TeV.pileup
        
Now everything should be setup for proper running, you can test it by doing the following

        rm test.root; xrd eoscms cat /eos/cms/store/cmst3/user/obondu/HH/signal/13TeV/MGraviton_900_HHtobbbb.hepmc.gz | gunzip | ./DelphesHepMC Snowmass_cards/delphes_card_Snowmass_NoPileUp.tcl test.root
        rm test.root; xrd eoscms cat /eos/cms/store/cmst3/user/obondu/HH/signal/13TeV/MGraviton_900_HHtobbbb.hepmc.gz | gunzip | ./DelphesHepMC Snowmass_cards/delphes_card_Snowmass_50PileUp.tcl test.root

### Run Delphes through lxbatch

A couple of utilities are here to run on lxbatch

## Setup the analysis

Yet to be written

## OLD INSTRUCTIONS BELOW

<ul>
  <li>run do_links (assume you have installed this repository as one subdirectory of Delphes-3.0.9/10)
    <ul>
      <li>./do_links.sh</li>
    </ul>
  </li>
  <li>setup (used to setup correct root version, compatible with .root creation on lxplus, and libraries directory while running)
    <ul>
      <li>source setup.sh</li>
      <li>source setup_afs.sh</li>
    </ul>
  </li>
  <li>compile (where BLA is the program defined in BLA.cc)
    <ul>
      <li>make BLA.exe</li>
    </ul>
  </li>
  <li>run
    <ul>
      <li>./BLA.exe</li>
      <li>./ntupleProducerVbfHHbbXX.exe -i ../data/delphes_output_vbfHH_MCHM4_0_06_13_GEN.root</li>
    </ul>
  </li>
</ul>






Delphes simulation:

    rm data/delphes_output_ggHHnew.root
    ./DelphesHepMC examples/delphes_card_CMS.tcl data/delphes_output_ggHHnew.root data/test-MR410_out.lhe.hepmc
    ./DelphesHepMC examples/delphes_card_CMS.tcl data/delphes_output_ttbar.root data/test-ttbar.lhe.hepmc
    ./DelphesHepMC examples/delphes_card_CMS.tcl data/delphes_output_vbfHHWWbb.root data/test-MR410_out.lhe.hepmc
    

Ntuple producer:

    cd LH-HH-Analysis
    source setup.sh 
    make ntupleProducerVbfHHbbXX.exe
    ./ntupleProducerVbfHHbbXX.exe -i ../data/delphes_output_ggHHnew.root
    ./ntupleProducerVbfHHbbXX.exe -i ../data/delphes_output_ttbar.root
    ./ntupleProducerVbfHHbbXX.exe -i ../data/delphes_output_vbfHHWWbb.root
