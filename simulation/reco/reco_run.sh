#! /bin/bash

#pathData=/home/mlf/zayunsna/data_storage/ext_gamma_thermalN_220916/jade
pathData=/home/mlf/zayunsna/data_storage/ext_gamma/ext_gamma_thermalN_220921/jade

subRun=$1

input=${pathData}/ext_gamma_thermalN_sub${subRun}.root

output=${pathData}/reco/Reco_ext_gamma_thermalN_sub${subRun}.root

echo " ---> jade_reco -i ${input} -o ${output} -m Q"
jade_reco -i ${input} -o ${output} -m Q

