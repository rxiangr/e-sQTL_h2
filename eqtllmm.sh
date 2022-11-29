#!/bin/bash
#SBATCH --job-name=makeGRM
#SBATCH --mail-type=FAIL                               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ruidong.xiang@agriculture.vic.gov.au      # Where to send mail
#SBATCH --ntasks=1                                      # Run a single task
#SBATCH --cpus-per-task=1                              # Run on a single or multiple CPUs
#SBATCH --nodes=1                                       # Run on a single Node
#SBATCH --mem=60gb                                     # Job memory request
#SBATCH --time=100:00:00                                # Time limit hrs:min:sec
#SBATCH --account=dbioanim1
#SBATCH --error=%x.%j.err
#SBATCH --output=%x.%j.out

#---variables

#---convertplink memery unit
let "memmb = $memg * 1000"

echo "Job of $pref started at"|paste -d ' ' - <(date)

#==========
cd $nodeDir
#==========
#---make temp file
mkdir gwatemp
if [ $? -gt 0 ]; then
       echo "ERROR; Failed to make temp files on comp node"
       exit 1
       fi
gwaout=gwatemp
#----loop through gene list (for those genes appeared in the pheno header file)
tail -n +2 $geneanno|awk -F '\t' -v chrn="$chrn" '$3==chrn{print$1}'|grep -f $phenopath/$phenheaderfn - > genelist
Ngene=`wc -l genelist|awk '{print$1}'`
for gene in `cat genelist`
  do ordergene=`grep -n $gene genelist|cut -d':' -f1`
     echo "Doing gene of $gene ($ordergene out $Ngene) at"|paste -d ' ' - <(date)
  #---copy pheno files to node/gene
  cp $phenopath/{$phenfn,$phenheaderfn} .
  if [ $? -gt 0 ]; then
       echo "ERROR; Failed copy files to comp node"
       exit 1
       fi

  #---copy cov files to node/gene
  cp $covpath/{$covfn,$qcovfn} .
  if [ $? -gt 0 ]; then
       echo "ERROR; Failed copy files to comp node"
       exit 1
       fi

 #---process phenodata based on fam
   trNum=`grep -n $gene $phenheaderfn|cut -d':' -f1`
   subset -i $phenfn -o pheno0.txt -s $plinkindir/$plinkfn.fam -c 1 -k 2
   reorder.scr pheno0.txt 1 $plinkindir/$plinkfn.fam 2|awk -v col="$trNum"  '{print $1,$2,$(col+2)}' > pheno.txt
   if [ $? -gt 0 ]; then
       echo "ERROR; Failed to sort pheno on comp node"
       exit 1
       fi

 #---process cov data
   subset -i $covfn -o cov0.txt -s $plinkindir/$plinkfn.fam -c 1 -k 2
   reorder.scr cov0.txt 1 $plinkindir/$plinkfn.fam 2 > cov.txt
   if [ $? -gt 0 ]; then
       echo "ERROR; Failed to sort cov on comp node"
       exit 1
       fi

  #---process qcov data
   subset -i $qcovfn -o qcov0.txt -s $plinkindir/$plinkfn.fam -c 1 -k 2
   reorder.scr qcov0.txt 1 $plinkindir/$plinkfn.fam 2 > qcov.txt
   cut -d ' ' -f1-15 qcov.txt > qcov.new.txt
   cut -d ' ' -f1-4,15 qcov.txt > qcov.new2.txt
   if [ $? -gt 0 ]; then
       echo "ERROR; Failed to sort qcov on comp node"
       exit 1
       fi
 
  #---gcta gwas
  /group/dairy/for_Ruidong_multitrait/rx01/bin/mta01102019/gcta_1.93.1beta/gcta64 --mlma --reml-no-constrain --reml-bendV  --reml-maxit 200 --bfile $plinkindir/$plinkfn --grm $GWgrmpath/${tiss} --covar cov.txt --qcovar qcov.txt  --pheno pheno.txt --out $gwaout/$gene --thread-num $corNum
    if [ $? -gt 0 ]; then
    echo "Failed attempt of gcta gwas of $gene with full cov"
    /group/dairy/for_Ruidong_multitrait/rx01/bin/mta01102019/gcta_1.93.1beta/gcta64 --mlma --reml-no-constrain --reml-bendV  --reml-maxit 200 --bfile $plinkindir/$plinkfn --grm $GWgrmpath/${tiss} --covar cov.txt --qcovar qcov.new.txt  --pheno pheno.txt --out $gwaout/$gene --thread-num $corNum   
      if [ $? -gt 0 ]; then
      echo "Failed attempt of gcta gwas of $gene reduced qcov1"
      /group/dairy/for_Ruidong_multitrait/rx01/bin/mta01102019/gcta_1.93.1beta/gcta64 --mlma --reml-no-constrain --reml-bendV  --reml-maxit 200 --bfile $plinkindir/$plinkfn --grm $GWgrmpath/${tiss} --covar cov.txt --qcovar qcov.new2.txt  --pheno pheno.txt --out $gwaout/$gene --thread-num $corNum
        if [ $? -gt 0 ]; then
        echo "Failed attempt of gcta gwas of $gene reduced qcov2"
        /group/dairy/for_Ruidong_multitrait/rx01/bin/mta01102019/gcta_1.93.1beta/gcta64 --mlma --reml-no-constrain --reml-bendV --reml-alg 1 --reml-maxit 200 --bfile $plinkindir/$plinkfn --grm $GWgrmpath/${tiss} --covar cov.txt --qcovar qcov.new2.txt  --pheno pheno.txt --out $gwaout/$gene --thread-num $corNum
          if [ $? -gt 0 ]; then
          echo "Failed attempt of gcta gwas of $gene reduced qcov2 and reml alg 1"
          /group/dairy/for_Ruidong_multitrait/rx01/bin/mta01102019/gcta_1.93.1beta/gcta64 --mlma --reml-no-constrain --reml-bendV --reml-alg 2 --reml-maxit 2000 --bfile $plinkindir/$plinkfn --grm $GWgrmpath/${tiss} --covar cov.txt --qcovar qcov.new2.txt  --pheno pheno.txt --out $gwaout/$gene --thread-num $corNum
           if [ $? -gt 0 ]; then
           echo "Failed attempt of gcta gwas of $gene reduced qcov2 and reml alg 2"
           fi 
          fi
        fi
      fi
    fi
  
     if [ $? -gt 0 ]; then
       echo "ERROR; Failed to copy gwas results of $gene to $gwaout"
       #exit 1
       fi
    echo "Job of $pref for $gene finished at"|paste -d ' ' - <(date)
done

