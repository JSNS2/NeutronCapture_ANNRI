#! /bin/bash

subRun=$1

jade=${PWD}
#input=/group/mlf/nu/data/HyoungKu_Jeon/ext_gamma_220914
input=/group/mlf/nu/data/HyoungKu_Jeon/ext_gamma/ext_gamma_thermalN_220921
output=${input}/jade

#thrs=${jade}/dat/thresholdsDataLS.root
HG=${jade}/dat/highGainFitsLS.root
#HG=${jade}/dat/highGainFitsLS_new2.root
ratio=${jade}/dat/gainRatiosLS.root


jade_process_data -i ${input}/ext_gamma_thermalN_${subRun}.root -o ${output}/ext_gamma_thermalN_sub${subRun}.root -g ${HG} -r ${ratio} 

