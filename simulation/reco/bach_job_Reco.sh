#!/bin/bash

pathHere=${PWD}
#pathData=/home/mlf/zayunsna/data_storage/ext_gamma_thermalN_220916/jade
pathData=/home/mlf/zayunsna/data_storage/ext_gamma/ext_gamma_thermalN_220921/jade

log=${pathData}/log

subNum=$1


for i in $(seq 0 ${subNum})
	do
		log=${pathData}/log/Log_gamma_sub${i}.log
		#log=${pathData}/log/Log_michelE_sub${sub}.log
		bsub -o ${log} -q s ./reco_run.sh ${i}
		#echo "./reco_run.sh ${runNum} ${i} -o ${log}"
		#./reco_run.sh ${i}
	done

