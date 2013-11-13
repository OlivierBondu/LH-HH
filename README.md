LH-HH
=====

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

Ntuple producer:

    cd LH-HH-Analysis
    source setup.sh 
    ./ntupleProducerVbfHHbbXX.exe -i ../data/delphes_output_ggHHnew.root