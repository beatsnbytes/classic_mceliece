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
    sed  -e "s+/$2+/$1+g" -e  "s/$2_$2/$1_$kernel_num/g" -e  "s+$2\*+$kernel_num\*+g" -e "s+$prev_mult\*+$prev_idx\*+g"  synd_kernel$2_$2.c > synd_kernel$1_$kernel_num.c
    echo "test_kernel_$1_$kernel_num"
done

# sed  -e "s+/$2+/$1+g" -e  "s/$2_1/$1_1/g" synd_kernel$2_1.c > synd_kernel$1_1.c
# echo "test_kernel_$1_1"

# start=$((1 + 1))
# prev_mult=$((4 - 1))
# for (( kernel_num=$start; kernel_num<=$1; kernel_num++ ))
# do
#     prev_idx=$(($kernel_num - 1)) 
#     sed  -e "s+/$2+/$1+g" -e  "s/$2_4/$1_$kernel_num/g" -e "s+$prev_mult\*+$prev_idx\*+g" -e  "s+4\*+$kernel_num\*+g"  synd_kernel$2_4.c > synd_kernel$1_$kernel_num.c
#     echo "test_kernel_$1_$kernel_num"
# done