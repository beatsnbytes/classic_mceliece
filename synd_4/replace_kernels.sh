#!/bin/bash

# Give as arg $1 the multiplicity you want to procude ans as $2 the multiplicity you currently have
for (( kernel_num=1; kernel_num<=$2; kernel_num++ ))
do
    sed  -e "s+/$2+/$1+g" -e  "s/$2_$kernel_num/$1_$kernel_num/g" synd_kernel$2_$kernel_num.c > synd_kernel$1_$kernel_num.c
    echo "test_kernel_$1_$kernel_num"
done

start=$(($2 + 1))
prev_mult=$(($2 - 1))
for (( kernel_num=$start; kernel_num<=$1; kernel_num++ ))
do
    prev_idx=$(($kernel_num - 1)) 
    gmem_idx=$(($kernel_num % 9)) 
    sed  -e "s+/$2+/$1+g" -e  "s/$2_$2/$1_$kernel_num/g" -e  "s+$2\*+$kernel_num\*+g" -e "s+$prev_mult\*+$prev_idx\*+g" -e"s/gmem1/gmem$gmem_idx/g"  synd_kernel$2_$2.c > synd_kernel$1_$kernel_num.c
    echo "test_kernel_$1_$kernel_num"
done