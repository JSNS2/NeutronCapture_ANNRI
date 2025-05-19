#!/bin/bash

verbose=0

for i in {1..1000..1}
do
	bsub -q s sh RunSingle.sh $i

done
