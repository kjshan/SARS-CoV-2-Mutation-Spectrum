#Download and combined data
#Nanoball data from  https://doi.org/10.17605/OSF.IO/8F6N9 
#UMI == Barcode in paper
###################################################################

#Mapping to genome

###################################################################

STAR  --runThreadN 30  --outFilterMultimapNmax 20 --genomeDir ./Vero/STAR_index/ --readFilesIn  all_1.fastq all_2.fastq   --outFileNamePrefix Vero_ --outFilterType BySJout  --alignSJoverhangMin 8 --outSJfilterOverhangMin 12 12 12 12 --outSJfilterCountUniqueMin 1 1 1 1 --outSJfilterCountTotalMin 1 1 1 1 --outSJfilterDistToOtherSJmin 0 0 0 0 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --scoreGapNoncan -4 --scoreGapATAC -4 --chimOutType WithinBAM HardClip --chimScoreJunctionNonGTAG 0 --alignSJstitchMismatchNmax -1 -1 -1 -1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 
STAR  --runThreadN 30  --outFilterMultimapNmax 20 --genomeDir ./SARS_CoV_2/STAR_index/  --readFilesIn  all_1.fastq all_2.fastq    --outFileNamePrefix SCV2_ --outFilterType BySJout  --alignSJoverhangMin 8 --outSJfilterOverhangMin 12 12 12 12 --outSJfilterCountUniqueMin 1 1 1 1 --outSJfilterCountTotalMin 1 1 1 1 --outSJfilterDistToOtherSJmin 0 0 0 0 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --scoreGapNoncan -4 --scoreGapATAC -4 --chimOutType WithinBAM HardClip --chimScoreJunctionNonGTAG 0 --alignSJstitchMismatchNmax -1 -1 -1 -1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 

###################################################################

#Keep reads uniquely mapped to SARS-CoV-2 genome

###################################################################

perl Virus_unique.pl Vero_Aligned.out.sam SCV2_Aligned.out.sam > CV_Unique.sam
samtools sort -@ 10 CV_Unique.sam -o CV_Unique.bam

###################################################################

#The read distribution in SARS-CoV-2 genome

###################################################################

nohup bedtools genomecov -d -split -ibam CV_Unique.bam > CV_Unique.depth &

###################################################################

#Call SNP (Backgroud)

###################################################################

mkdir All_mpileup
samtools mpileup --max-depth 0 --output-BP --output-QNAME -a --reference /Dell/Dell9/shankj/Cov/RTmutation/SARS-CoV-2/SARS_CoV_2/NC_045512.2.fna CV_Unique.bam -o CV_Unique.mpileup 
cd All_mpileup
split -l 30 ../CV_Unique.mpileup -d -a 4
rm  ../CV_Unique.mpileup
ls|grep -v list > list.txt
parallel  -j 10 "perl mpileup.pl  ./{} > ./{}.out" :::: list.txt
cat ./*.out > ../SCV2.Unique.summarize
cd ../
perl total_position_mismatch_rate.pl SCV2.Unique.summarize > Guassi.txt

###################################################################

#Junction reads

###################################################################

perl Junction_read.pl CV_Unique.sam
samtools sort -@ 10  CV_Unique.sam.Junction.read -o CV_Unique.Junction.read.bam

###################################################################

#Each barcode contains how many read pairs

###################################################################

perl Junction_read_bed.pl CV_Unique.sam.Junction.read > CV_Unique.Junction.read.bed 
perl Count_Junction.pl CV_Unique.Junction.read.bed > CV_Unique.Junction.read.bed.count
perl SCV2_unique_UMI.pl CV_Unique.sam.Junction.read CV_Unique.Junction.read.bed CV_Unique.Junction.read.bed.count > SCV2_Unique_UMI.txt

###################################################################

#Mismatches in junction reads

###################################################################

mkdir JCR_mpileup
samtools mpileup --max-depth  0 --output-BP --output-QNAME -a --reference /Dell/Dell9/shankj/Cov/RTmutation/SARS-CoV-2/SARS_CoV_2/NC_045512.2.fna CV_Unique.Junction.read.bam -o CV_Unique.Junction.read.mpileup 
cd JCR_mpileup
split -l 30 ../CV_Unique.Junction.read.mpileup -d -a 4
rm ../CV_Unique.Junction.read.mpileup
cd Sam_split
perl ../Sam_split.pl ../CV_Unique.Junction.read.sam ../SCV2_Unique_UMI.txt > list.txt
parallel -j 50 "samtools sort {} -o ./{}.bam" :::: list.txt
parallel -j 50 "samtools mpileup -d 0 --reference /Dell/Dell9/shankj/Cov/RTmutation/SARS-CoV-2/SARS_CoV_2/NC_045512.2.fna --output-QNAME --output-BP {}.bam -o ./{}.mpileup.default" :::: list.txt
parallel -j 50 "perl ../Mutation_Rate.pl {}.mpileup.default " :::: list.txt > ../Consensus_read.out #(ATCG content in consensus sites, 26257-26283 should be discarded)

perl Coverage_consensus_reads.pl ../Consensus_read.out > ../Coverage_consensus_reads.txt

cd ../
perl Mutation_count.pl SCV2.Unique.summarize SCV2.Unique.Perl.JCR.summarise > Vero_Total_Barcode_SNP.txt

###################################################################

#The positions in SARS-CoV-2 genome, covered by junction reads

###################################################################

parallel  -j 10 "perl Potential_mutation_site.pl ../SCV2_Unique_UMI.txt {} " :::: list.txt > ../Potential_mutation_site.txt 

###################################################################

#Consenus VS Unconsensus

###################################################################
samtools mpileup --max-depth  0 --output-BP --output-QNAME -a -Q 0 -B --adjust-MQ 0 --reference /Dell/Dell9/shankj/Cov/RTmutation/SARS-CoV-2/SARS_CoV_2/NC_045512.2.fna CV_Unique.Junction.read.bam -o CV_Unique.Junction.read.mpileup 
cd JCR_mpileup
split -l 30 ../CV_Unique.Junction.read.mpileup -d -a 4
rm ../CV_Unique.Junction.read.mpileup
ls |grep -v list > list.txt
parallel  -j 10 "perl Unconsensus.pl   ./{} " :::: list.txt > ../Unconsensus 
parallel  -j 10 "perl Consensus.pl     ./{} " :::: list.txt > ../consensus 
cd ../
perl Consensus2.pl consensus > Consensus.result
perl Unconsensus2.pl Unconsensus > Unconsensus.result


###################################################################

#Small Indel

###################################################################

perl UMI_indel.pl CV_Unique.sam.Junction.read > UMI_indel.sam
samtools view -@ 30 -bT /Dell/Dell13/shankj/projects/Cov/NC_045512.2.fna UMI_indel.sam > UMI_indel.bam
samtools sort -@ 30 UMI_indel.bam -o UMI_indel.sort.bam
rm UMI_indel.bam

#UMI
perl UMI.pl UMI_indel.sam > JCR_indel.txt2
perl UMI2.pl JCR_indel.txt2 > JCR_indel.txt2.gt1

#split Sam by UMI
mkdir Sam_Split
cd Sam_Split
perl Sam_Split.pl ../UMI_indel.sam2 JCR_indel.txt2.gt1
ls *.sam|grep -v list > list.txt

#call indel
parallel  -j 50 "samtools sort  {} -o ./{}.bam" :::: list.txt
parallel  -j 50 "samtools mpileup -d 0  --reference ./NC_045512.2.fna  --output-QNAME  --output-BP  {}.bam -o  ./{}.mpileup.default" :::: list.txt
parallel  -j 50 "perl Indel_information.pl {}.mpileup.default " :::: list.txt > ../UMI_indel.information

cd ../
perl UMI_indel_information_PCR.pl  > UMI_indel.information.PCR.duplication

awk '{if ($13>=2 && $3==$14) print}' UMI_indel.information.PCR.duplication > Consensuse.Indel.txt
    1                   2                    3                         4                    5      6      7        8        9       10          11      12          13               14
#"Barcode","Barcode_All_Read_pair_number","Covered_reads_number","Covered_readpair_number","Pos","Indel","Alt","Indel_len","Dis","Reads_Pos","Reads","Reads_PCR","non_PCR_number","SNP_read_number"

#Background
awk '{if($6~/M/ && ( $6~/I/ || $6~/D/)  && $6!~/N/  && $12 eq "NH:i:1") print }' CV_Unique.sam > JCR_indel.background.sam

samtools view -@ 30 -bT /Dell/Dell13/shankj/projects/Cov/NC_045512.2.fna JCR_indel.background.sam > JCR_indel.background.bam
samtools sort -@ 30 JCR_indel.background.bam -o JCR_indel.background.sort.bam
rm JCR_indel.background.bam
samtools mpileup -d 0  --reference /Dell/Dell13/shankj/projects/Cov/NC_045512.2.fna  --output-QNAME  --output-BP JCR_indel.background.sort.bam -o JCR_indel.background.mpileup.default
perl Indel_information.pl JCR_indel.background.mpileup.default > ./UMI_indel.information.background.default

          1                   2                    3       4      5      6         7        8        9       10          11      12          13               14
#"Covered_reads_number","Covered_readpair_number","Pos","Indel","Alt","Indel_len","Dis","Reads_Pos","Reads","Reads_PCR","non_PCR_number","SNP_read_number"


perl Polymorphism_Indel.pl > Polymorphism_Consensuse.Indel.txt
