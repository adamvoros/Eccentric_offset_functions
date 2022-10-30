#!/bin/sh

OPTSN=./optsn

for k in `seq -0.4 0.2 0.5` ; do
    for h in `seq -0.4 0.2 0.5` ; do
	A=($($OPTSN -e $k,$h -1 -n 4 | sed -e 's/#.*//'))
	printf "%6.1f %6.1f %6.4f %6.4f %6.4f %6.4f\n" $k $h ${A[*]} 
    done
done

