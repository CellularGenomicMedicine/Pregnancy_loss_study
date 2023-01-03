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
log_file=$3
script=$4
process=$5
analyses=$6
gtypemodulator_window=$7
family=$8
SibPattern=$9

cd $logs/error-logs

date=$(date '+%F_%H%M')
if [ ! -f $PID.finished.txt ]; then
 if [ ! -e $PID.BUSY ]; then
  touch $PID.BUSY
  if [ -e $PID.BUSY ]; then
   echo $date 'Start' $process $SibPattern >> $logs/$log_file.log;
  
   sh $script/analyses/$process/$process.sh $PID $logs $script $process $analyses $gtypemodulator_window $family $SibPattern
  
   SRC=$?
   if [ $SRC != 0 ]; then  
    echo $date "Error ($SRC) in" $process $SibPattern >> $logs/$log_file.log; 
   fi
  else
   echo $date 'Could not create PID' $process $SibPattern >> $logs/$log_file.log;
   exit 1 
  fi
 else
  echo $date 'Already created PID' $process $SibPattern >> $logs/$log_file.log;
  exit 1 
 fi
else
 echo $date 'Already finished' $process $SibPattern >> $logs/$log_file.log;
fi

exit
