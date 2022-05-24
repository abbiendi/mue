# MuE
C++ interface to MESMER NLO MC generator (https://github.com/cm-cc/mesmer)
MUonE Fast Detector Simulation and Analysis

There are three subdirectories:

- “code” contains all the code and the C++ data format definition for the MESMER MC generator; 

- “test” contains example scripts:
   -- test_mesmer_MuE.sh  generate MESMER events and does fast simulation and analysis producing output histograms and/or trees
   -- test_MuE.sh         does fast simulation and analysis of pre-generated MESMER events
 
- “writer” contains code to run the MESMER generator and write out the MC events with defined ROOT data format

The examples use the 'embedded' mode for event generation with MESMER, which allows to run the (FORTRAN) generator and the (C++) fast simulation and analysis in one step. 
Previous versions of the code used two separate processes for the generation and the fast simulation and analysis and are not supported anymore, although they could still be working.

This FastSimulation and Analysis code has been used for the analysis contained in the LoI and studies of systematic errors (see details in arXiv:2201.13177)

The structure of the MESMER MC events is described by the class Event in the header https://gitlab.cern.ch/muesli/nlo-mc/mue/-/blob/master/code/MuEtree.h . It consists of:

    Int_t RunNr;
    Long_t EventNr;
    Double_t wgt_full, wgt_norun, wgt_lep, wgt_LO, wgt_NLO;
    static const uint NC = 11;
    std::array<Double_t, NC> coef;
    Particle muin;
    std::vector<Particle> fspart;

    RunNr    : run number
    EventNr  : event number
    wgt_full : main event weight
    wgt_norun: weight without running of alpha
    wgt_lep  : weight with only leptonic running of alpha
    wgt_LO   : Leading Order weight
    wgt_NLO  : Next-to-Leading Order weight
    coef     : array of 11 coefficients needed for the reweighting
    muin     : particle representing the incoming muon
    fspart   : std::vector containing the final state particles

    particles are defined by the class Particle which contains:

    Short_t pdgId;
    Double_t px;
    Double_t py;
    Double_t pz;

    from these quantities the energy and the mass can be returned by functions:
    Double_t E() const;
    Double_t M() const;

In the same header MuEtree.h the class Setup contains the settings of all the generator parameters and counters of the generated statistics (number of generated weights, their sums and sums of the squares, cross section, etc…).
In particular the class MCPara contains all the generator settings, the class MCstat all the counters and quantities calculated accumulating the events.

The main program https://gitlab.cern.ch/muesli/nlo-mc/mue/-/blob/master/code/mesmer_MuE.cc is the steering code which does event generation with MESMER, embedding its function calls, and at the same time the fast simulation and analysis. 
In the subdirectory test, the script test_mesmer_MuE.sh compiles everything and runs a test job, by generating and processing 100K events.
The test job creates a working directory producing the histogram file called results.root and optionally a tree called outtree.root, together with few test plots of kinematic distributions and resolutions.

The main program https://gitlab.cern.ch/muesli/nlo-mc/mue/-/blob/master/code/MuE.cc is used to process pre-generated MESMER events, allowing for multiple input files. It reads in multiple files and the related infos in the Setup, checking their consistency, and loops over the events doing the fast simulation and analysis, as in the previous example.

MESMER events can be produced and saved in output ROOT format without any other processing by the code 
https://gitlab.cern.ch/muesli/nlo-mc/mue/-/blob/master/writer/mesmer_writer.cc, which is used in the example test_mesmer_writer.sh.

All the examples check for the existence of the MESMER installation directory and in case download the code from github and compile it.

My suggestion to start with is trying to run the test example and play with it.
In short: on lxplus.cern.ch

```
$ git clone https://gitlab.cern.ch/muesli/nlo-mc/mue.git
$ cd mue/test
$ ./test_mesmer_MuE.sh
```


