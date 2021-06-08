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
perl SingleToMultiple.pl
mkdir Astr
cd Astr
perl /Dell/Dell13/shankj/bin/FastML.v3.11/www/fastml/FastML_Wrapper.pl --MSA_File ../result.txt --seqType  NUC --outDir ../Astr > fastML.out 
grep  N1$ -A 1  seq.joint.txt > ../Astr.fas
cat ../Astr.fas ../align.fas > ../Astr_Align.fas
cd ..

#call SNP, ancestral sequence as reference genome 
perl MSA_SNP2.pl ./Astr_Align.fas N1 > ${file%.fas}_Astr_Align.SNP
rm result.txt
