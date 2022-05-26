#!/bin/sh

wkdir=`pwd`
gendir=${wkdir}/../mesmer

if [ ! -d ${gendir} ] 
then
    mkdir ${gendir}
    cd ${gendir}
    echo "Installing MESMER at directory: " ${gendir} 
    git clone https://github.com/cm-cc/mesmer.git -b v1.1.0 
    make
else
    cd ${gendir}
    echo "MESMER installation directory is: " ${gendir}
    make
fi

cd ${wkdir}
