#!/bin/sh

writer=`pwd`
gen=${writer}/../mesmer
src=${writer}/../code
test=${writer}/../test

# MESMER setup
source ${test}/mesmer_setup.sh

MAIN=mesmer_writer

ROOTINCDIR=`$ROOTSYS/bin/root-config --incdir`
ROOTLIBS=`$ROOTSYS/bin/root-config --cflags --libs`

# working directory
JOBDIR=test_mesmer_writer_`date +"%y-%m-%d_%T"`
echo "Working directory is: " ${JOBDIR}
mkdir ${JOBDIR}
cd ${JOBDIR}

echo "Compiling ..."

rootcling -f MuEtreeWDict.C -I${src} MuEtree.h Utils.h -I${writer} MuEtreeWLinkDef.h

OPTIONS="-O2 -Wall"
#OPTIONS="-g -Wall"

g++ ${OPTIONS} -shared -fPIC -o libMuEtree.so -I${src} MuEtreeWDict.C ${src}/MuEtree.cc ${src}/Utils.cc ${ROOTLIBS}

g++ ${OPTIONS} -o ${MAIN}.exe -I${src} -I${ROOTINCDIR} ${writer}/${MAIN}.cc MuEtreeWDict.C ${src}/MuEtree.cc ${src}/Utils.cc ${ROOTLIBS} -L${gen}/ -lmesmerfull -lX11 -lgfortran

cat > mesmer.cards <<!
mode weighted
nev 100000
nwarmup 10000
seed 42
sync 1
ord alpha
Ebeam 160.d0
bspr  3.75d0
extmubeam no
Qmu 1
Eemin 0.2d0
thmumin 0.
thmumax 100.
themax 100.
path test-run
run
!

echo "Running " ${MAIN}.exe
time ./${MAIN}.exe  mesmer.cards  mesmer_events.root
#time ./${MAIN}.exe  mesmer.cards  mesmer_events.root  >&  ${MAIN}.log 2>&1

rm *.exe *.C

echo "Done."
