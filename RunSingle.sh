#!/bin/bash

source /home/mlf/cdshin/JADE/MC/setup.sh

num=$1

verbose=0

/home/mlf/cdshin/RAT/jsns2/externalVer/ext_gamma_gen/JSNS2_ANNRI_Gd/obj/reGeneration_nGd_gamma $num test__ $verbose

