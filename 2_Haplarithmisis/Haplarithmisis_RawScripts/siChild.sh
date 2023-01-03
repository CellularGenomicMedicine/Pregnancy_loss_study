#!/bin/bash
# kasper.derks
#$ -S /bin/bash
#$ -cwd
. /etc/profile
#----------------------------------------------------------------------
source /etc/profile.d/modules.sh
module load gcc/4.8.1
module load ngs_prd
module load R/3.3.3
module load samtools/1.3

PID=$1
logs=$2
script=$3
process=$4
analyses=$5
gtypemodulator_window=$6
family=$7
SibPattern=$8

Rscript $script/analyses/$process/$process.R $PID $logs $script $analyses $gtypemodulator_window $family $SibPattern

