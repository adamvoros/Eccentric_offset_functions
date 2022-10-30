#!/bin/sh

# This script shows an example how 5 phases from the total of n=8 can be 
# derived involving the ./optsn code when other 3 of the phases are kept
# fixed. In mode `-1`, the ./optsn utilizes a downhill simplex minimization
# for the phases that are not fixed to minimize the phase volume U in the
# (k,h) parameter space. These phases are initiated for random positions
# for 1000 times (for this particular case, 1000 minimizations are sufficient
# to find the global minimum). The output is then sorted by the s/n ratio (last
# column of the output, following the '#' mark), thus the last few lines
# printed by the `sort` command will show the optimal phases. These last
# lines should be something like:
# 0.12645 0.20000 0.30000 0.40000 0.42227 0.56910 0.56911 0.85390 #  123.418

OPTSN=./optsn

n=8

for s in `seq 1000` ; do 
	$OPTSN -e 0,0 -1 -n $n -s $s -p 0.2,0.3,0.4 
done | \
sort -k $((n+2)),$((n+2))g



