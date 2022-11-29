#!/bin/bash -l
#SBATCH --tasks-per-node=1 # //PBS number of nodes
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G # //PBS memory
#SBATCH --time=72:00:00 # //PBS wall time
#SBATCH --job-name=byr3 # //PBS job name
#SBATCH --account=dbioanim1  # //PBS group
#SBATCH -o job%j.o # //creats PBS *o* file, '%j' is the job id
#SBATCH -e job%j.e # //creates PBS *e* file, '%j' is the job id
#SBATCH --mail-user=ruidong.xiang@ecodev.vic.gov.au # //PBS email account
#SBATCH --mail-type=FAIL # //PBS when to send an email: BEGIN, END, FAIL, REQUEUE or ALL

#---variables

#run BayesR:
module purge
module load intel/2021b
echo $byrDirbayesR3 $pref.param.txt  >> $logFile
$byrDir/bayesR3 $pref.param.txt  2>&1 >> $logFile

if [ $? -gt 0 ]; then
    echo "bayesR3 failed"
    echo "Output files retained on $nodeDir"
    exit 1
fi

#${binDir}zipOnNodeAndCopyBack.scr $nodeDir $outDir *_aveEffects.csv
#${binDir}zipOnNodeAndCopyBack.scr $nodeDir $outDir *_Model.chn*.csv
cp *para* $outDir
cp *map.gz $outDir
cp *.log $outDir
if [ $? -gt 0 ]; then
    echo "Failed to copy output files"
    exit 1
fi

echo "byr run fihished and copied output files to" $outDir

cd ../
rm -r $nodeDir

