#!/bin/bash

subNum=$1



for i in $(seq 0 ${subNum})
	do
		bsub -o test.log -q s ./proc_run_normal.sh ${i} 
	done

