##################################################################################################################################################

#Call polymorphisms
#SARS-CoV-2,H1N1 and MERS from patients

##################################################################################################################################################
file=$1 #RSV_B.fas 
refname=$2 #KF893260.fas

#Change fasta format
bioawk -c fastx '{print ">"$name"\n"$seq}' ${file%.fas}.fasta > $file
bioawk -c fastx '{print ">"$name"\n"$seq}' ${refname%.fas}.fasta > $refname

mkdir Seq_Split 
mkdir muscle 

#Every sequence align to 
cd ./Seq_Split 
perl seq_split.pl ../$file ../$refname 
ls *.fas |grep -v "list" > ./list.txt 
parallel  -j 5 "muscle -in ./{} -out ../muscle/{}.out" :::: ./list.txt 

#combine alignment results
cd ../muscle/ 
ls *.fas.out |grep -v "list"  > list.txt 
parallel  -j 60 "perl align_combine.pl ./{}  ${refname%.fas} > ./{}.align" :::: ./list.txt  
cat *.align > ../total.align 
cat ../$refname ../total.align  > ../align.fas 
rm  ../total.align 

#call SNP, the genome you provided  as reference genome 
perl MSA_SNP2.pl ../align.fas ${refname%.fas} > ../${refname%.fas}.SNP.txt 

#construct ancestral sequence using FastML
cd ../ 
perl /Dell/Dell13/shankj/projects/Cov/Genome/NegativeStrand/Kejia/HKU/SingleToMultiple.pl ./SARS_nopangolin.muscle.out

mkdir Astr
cd Astr
perl /Dell/Dell13/shankj/bin/FastML.v3.11/www/fastml/FastML_Wrapper.pl --MSA_File ../SARS-CoV_10.out --seqType  NUC --outDir ../Astr > fastML.out  &
grep  N1$ -A 1  seq.joint.txt > ../Astr.fas
cat ../Astr.fas ../align.fas > ../Astr_Align.fas

##################################################################################################################################################

#Call polymorphisms in the Fig6-7
#SARS-CoV-2,SARS and MERS from patients

##################################################################################################################################################

cd /Dell/Dell13/shankj/projects/Cov/Genome/Fig6/20210710
mkdir SCV2
cd SCV2
muscle -in SCV2.fasta -out SCV2.muscle.out
perl SingleToMultiple.pl ./SCV2.muscle.out
mkdir Astr
cd Astr
perl /Dell/Dell13/shankj/bin/FastML.v3.11/www/fastml/FastML_Wrapper.pl --MSA_File ../result.fas --seqType  NUC --outDir ../Astr > fastML.out 
perl /Dell/Dell13/shankj/projects/Cov/Genome/Fig6/MERS/Astr/Astr_Branch.pl SCV2 > SCV2_SNP.count


cd /Dell/Dell13/shankj/projects/Cov/Genome/Fig6/20210711/SARS/Astr/
perl /Dell/Dell13/shankj/projects/Cov/Genome/Fig6/MERS/Astr/Astr_Branch.pl SARS > ../SARS_SNP.count
cd /Dell/Dell13/shankj/projects/Cov/Genome/Fig6/20210709/SCV2/Astr/
perl /Dell/Dell13/shankj/projects/Cov/Genome/Fig6/MERS/Astr/Astr_Branch.pl SCV2 > /Dell/Dell13/shankj/projects/Cov/Genome/Fig6/20210709/SNP/SCV2_SNP.count
cd /Dell/Dell13/shankj/projects/Cov/Genome/Fig6/20210709/MERS/0710/Astr/
perl /Dell/Dell13/shankj/projects/Cov/Genome/Fig6/MERS/Astr/Astr_Branch.pl MERS > /Dell/Dell13/shankj/projects/Cov/Genome/Fig6/20210709/SNP/MERS_SNP.count
cd /Dell/Dell13/shankj/projects/Cov/Genome/Fig6/20210709/SNP/
cat *.count|grep -v "ranch" > SNP
