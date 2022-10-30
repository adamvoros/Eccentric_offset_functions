#!/bin/sh

# This script calculates the phases for N_obs=5, 6 or 7.

OPTSN=./optsn

for n in 5 6 7 ; do 

	for s in `seq 100` ; do 
		$OPTSN -e 0,0 -1 -n $n -s $s 
	done | \
	sort -k $((n+2)),$((n+2))g | \
	tail -n 1 

done



