#!/bin/sh

test=`pwd`
gen=${test}/../mesmer
src=${test}/../code

# MESMER setup
source ./mesmer_setup.sh

MAIN=MuE

ROOTINCDIR=`$ROOTSYS/bin/root-config --incdir`
ROOTLIBS=`$ROOTSYS/bin/root-config --cflags --libs`

# working directory
JOBDIR=test_MuE_`date +"%y-%m-%d_%T"`
echo "Working directory is: " ${JOBDIR}
mkdir ${JOBDIR}
cd ${JOBDIR}

echo "Compiling ..."
rootcling -f MuEtreeDict.C -I${src} MuEtree.h MuEana.h MuEtreeLinkDef.h

OPTIONS="-O2 -Wall"
#OPTIONS="-g -Wall"

g++ ${OPTIONS} -I${src} -I${ROOTINCDIR} ${src}/${MAIN}.cc MuEtreeDict.C ${src}/MuEtree.cc ${src}/Inputs.cc ${src}/Analysis.cc ${src}/FastSim.cc ${src}/Utils.cc ${src}/ElasticState.cc ${src}/Histos.cc ${src}/dalpha.cc ${ROOTLIBS} -L${gen} -lmesmerfull -lgfortran -lX11 -o ${MAIN}.exe

#exit

################################################################################
### FASTSIM input configuration ################################################
################################################################################
cat > fastSim.cfi <<!
<cfi>
0        # detector resolution model: 0=simplest-2par; 1=Antonio's 3 par
0        # 0:default; 1:MS only on X plane; 2: polar angle smearing
0        # true=twosteps: apply multiple scattering in model 0 and intrinsic resolution in a second step
0.0425   # total material thickness for model 0 (in X0) / 1.5cm Be = 0.0425 X0
0.02     # intr.ang.resol. for model 0 (in mrad)
0.       # signed bias of station length (um) 
0.       # uncertainty of station length (um)
</cfi>
!
################################################################################

################################################################################
### ANALYSIS input configuration ################################################
################################################################################
cat > analysis.cfi <<!
<cfi>
0       # bool makeTree; 1(0) = do(not) produce the output Tree
1       # bool doTemplates; 1(0) = do (not) produce 2D template histos
50.     # max Theta for event selection
32.     # max Theta for template histograms
0       # 1(0): do(not) make angle correlation plots in fine bins
2       # parameterization of hadronic running: 0:pol2; 1:LL; 2:LLmod
0       # 0:nominal Lumi; 1:LowLumi; 
5       # range around expected average (as number of sigmas on each axis)
2       # number of divisions within one sigma interval
0       # bool doEnergyScale; 1(0) = do (not) produce plots for E scale calibration
</cfi>
!
################################################################################

################################################################################
### MAIN input configuration ################################################
################################################################################
cat > input.cfi <<!
<cfi>
${test}/mesmer_input.txt # string input_dirs_file; // file containing a list of directory paths where to look for NLO MC events
100000          # long n_events; // events to be processed (N>0:N; 0:all; <0: read histograms from existing results)
mesmer_events   # string input file name
results.root    # string histo_ifname; // path to input histo file (kinematical distributions when n_events<0)
MuE             # string output_dir; // output directory name
fastSim.cfi     # string fastSim_ifname; // FastSim cfg file
analysis.cfi    # string analysis_ifname; // Analysis cfg file
</cfi>
!

echo "Running ..."

time ./${MAIN}.exe  input.cfi
#time ./${MAIN}.exe  input.cfi  > ${MAIN}.log 2>&1

rm *.exe *.C
echo "Done."
