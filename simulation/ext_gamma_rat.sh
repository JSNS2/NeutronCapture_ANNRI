#!/bin/bash


for i in {0..1000..1}
do
	#sub=$(printf "%05d\n" $i)
	sub=$(printf "%d\n" $i)
	cp ext_gamma.mac ./mac/ext_gamma.mac
		`sed s/NumT/${sub}/ ./mac/ext_gamma.mac > ./mac/ext_gamma_$i.mac`
	echo "${i}th macro was created "
	cp run.sh ./script/run_$i.sh
	echo "${i}th Script file was created "
	bsub -q l ./script/run_$i.sh $i
done


