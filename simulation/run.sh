#!/bin/bash

#MAINDIR=/home/mlf/zayunsna/official_work/generator/JSNS2_eventGenerator_pack/thermal_n/externalVer
MAINDIR=${PWD}
MACDIR=${MAINDIR}/mac
#OUTDIR=/group/mlf/nu/data/HyoungKu_Jeon/thermal_neutron
#OUTDIR=/group/mlf/nu/data/HyoungKu_Jeon/thermal_neutron_220816
#OUTDIR=/group/mlf/nu/data/HyoungKu_Jeon/thermal_neutron_220831
OUTDIR=/group/mlf/nu/data/HyoungKu_Jeon/ext_gamma/ext_gamma_thermalN_220921


#OUTDIR=/group/mlf/nu/data/HyoungKu_Jeon/Birks_test/thermalN/Original

SEQ=$1

rat $MACDIR/ext_gamma_${SEQ}.mac -o $OUTDIR/ext_gamma_thermalN_${SEQ}.root

