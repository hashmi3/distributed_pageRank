#!/bin/bash

rm slurm-*
#make triangle_counting
#make triangle_counting_parallel
make USE_INT=1 page_rank_parallel
#make USE_INT=1 page_rank
#make page_rank_parallel
sbatch submit.sh
squeue
echo "################################"
#sleep 2 && less sl*

