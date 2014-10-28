#!/bin/sh
. ~whiter/.bashrc

# SGE OPTIONS now in qsub_options.txt

#funcon_proc.py args
project=$1
atlas=$2
suffix=$3
subblock=$4
smooth_opt=$5
filter_opt=$6
scrub_opt=$7
corr_opt=$8
subidx=$9

/home/despo/whiter/gitrepos/podNtools/funccon_code/funccon_proc.py $project $atlas $suffix $subblock $smooth_opt $filter_opt $scrub_opt $corr_opt $subidx
