#!/bin/bash

rm kat

for opt_level in 0 1 2 3
do
	make kat ARCH=arm OPT=$opt_level VEC=y
	mv kat timing_binaries/timing_arm/vec/kat_O$opt_level

	make kat ARCH=arm OPT=$opt_level VEC=n
        mv kat timing_binaries/timing_arm/no_vec/kat_O$opt_level

	make kat OPT=$opt_level VEC=y
        mv kat timing_binaries/timing_x86/vec/kat_O$opt_level

        make kat OPT=$opt_level VEC=n
        mv kat timing_binaries/timing_x86/no_vec/kat_O$opt_level


done


