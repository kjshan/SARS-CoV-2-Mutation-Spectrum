WORKDIR=$1    #H2O2_20min-2; sample name==work dir name
REFIDX=$2     #/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ours/H2O2_Test/Clean/STAR/Yeast_virus_STAR_index/
GENOME=$3     #/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ours/H2O2_Test/Clean/STAR/Yeast_virus_Genome.fasta
SCRIPTDIR=$4  #/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ours/H2O2_Test/Clean/STAR/data/bin/
MRNA=$5       #/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/Reference/mRNA.bed
THREAD=$6     #THREAD
#########################################################################################################################

#########################################################################################################################

#Map to Genome
#100% Match reads

#########################################################################################################################

cd ${WORKDIR}
#STAR  --runThreadN ${THREAD} --genomeDir ${REFIDX} --readFilesIn  ${WORKDIR}_1.fq.gz ${WORKDIR}_2.fq.gz --outFileNamePrefix ${WORKDIR} --outFilterType BySJout --readFilesCommand zcat
#perl ${SCRIPTDIR}/fragment_barcode.pl ${WORKDIR}Aligned.out.sam > ${WORKDIR}.noIndel.150M.sam

#########################################################################################################################

#Fragment Barcode and call SNP by samtools

#########################################################################################################################

#perl ${SCRIPTDIR}/Readpair_Barcode.pl ./${WORKDIR}.noIndel.150M.combined.sam > list.txt
#samtools sort -@ ${THREAD} ${WORKDIR}.noIndel.150M.combined.sam.sam2 -o ${WORKDIR}.noIndel.150M.combined.bam
#rm ${WORKDIR}.noIndel.150M.combined.sam.sam2 &
#samtools mpileup -d 0  --reference ${GENOME}  --output-QNAME  --output-BP  ${WORKDIR}.noIndel.150M.combined.bam -o  ./${WORKDIR}.noIndel.150M.combined.mpileup

#########################################################################################################################

#SNP Extraction
#split mpileup

#########################################################################################################################

#mkdir mpileup
#cd mpileup
#split -l 100000 ../${WORKDIR}.noIndel.150M.combined.mpileup -d -a 8 
#ls |grep -v 'list'  > list.txt
#parallel  -j ${THREAD} "perl ${SCRIPTDIR}/Readpair_Barcode_Mpileup.pl ../list.txt {}" :::: list.txt >> ../Coverage.default
#cat *.SNP > ../${WORKDIR}.Consensuse.Mismatch
##         1       2          3             4             5          6
##"Chr_Position\tRef > Alt\tAlt_number\tAlt_fraction\tRef_fraction\tCover\n";
#rm ./${WORKDIR}.noIndel.150M.combined.mpileup &

#########################################################################################################################

#Polarize the mRNA errors

#########################################################################################################################

#cd ../
#perl ${SCRIPTDIR}/Barcode_SNP_2.pl ./${WORKDIR}.Consensuse.Mismatch ${MRNA} > ./${WORKDIR}.Consensuse.Mismatch.Strand
perl ${SCRIPTDIR}/mRNA_Cover.pl ./Coverage.default ${MRNA} > ./Coverage.default.mRNA
perl ${SCRIPTDIR}/SNP_Strand.pl ./Coverage.default.mRNA ./${WORKDIR}.Consensuse.Mismatch.Strand > ./${WORKDIR}.Strand.SNP

cp ./${WORKDIR}.Strand.SNP /Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ours/SNP/bin2/ #Which dir you want to copy to
#########################################################################################################################

#Somatic minus
#Just some examples, you should change your files to obtain your results

#########################################################################################################################

#perl /Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ours/bin/Somatic_mutation.pl ./Ctrl2.Consensuse.Mismatch.Strand ./H2O2-1h-2.Consensuse.Mismatch.Strand > R2.somatic
#perl /Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ours/bin/Somatic_mutation_minus.pl R2.somatic ./Ctrl2.Strand.SN ./H2O2-1h-2.Strand.SNP



#########################################################################################################################

#Consensus VS Unconsensus

#########################################################################################################################

#samtools mpileup -d 0  -Q 0 -B --adjust-MQ 0  --reference ${GENOME}  --output-QNAME  --output-BP  ${WORKDIR}.noIndel.150M.combined.bam -o  ./${WORKDIR}.noIndel.150M.combined.mpileup2
#mkdir ConsensusVS
#cd ./ConsensusVS/
#split -l 100000 ../${WORKDIR}.noIndel.150M.combined.mpileup2 -d -a 8 
#ls |grep -v 'list'  > list.txt
#parallel  -j 15 "perl ${SCRIPTDIR}/ConsensusVS.pl ../list.txt {}" :::: list.txt 
#cat *.Consensus > ../Consensus.txt
#cat *.Unconsensus > ../Unconsensus.txt
#perl ${SCRIPTDIR}/ConsensusVS2.pl ../Consensus.txt > ../Consensus_mismatch.txt
#perl ${SCRIPTDIR}/ConsensusVS2.pl ../Unonsensus.txt > ../Consensus_mismatch.txt
#rm ../${WORKDIR}.noIndel.150M.combined.mpileup2


#########################################################################################################################

#PS: We called mRNA errors in yeast CirSeq by the scripts from https://github.com/LynchLab/TranscriptErrors
