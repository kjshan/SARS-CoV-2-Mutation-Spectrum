#######################################

#mRNA errors in HeLa cell CirSeq

#######################################

#The Scripts of CirSeq from https://andino.ucsf.edu/toolsandprotocols
#So, we just use the sam files to call mRNA errors.
#We use "samtools sort" to sort sam files to bam files
#You should change ${1} to your file name

samtools mpileup  ${1}.bam  -Q 20 --max-depth  0  --reference Homo_sapiens.GRCh37.dna.primary_assembly.fa -o ${1}.mpileup
perl mpileup_out.pl ${1}.mpileup > ./${1}.mpileup.out 
perl mRNA_error.pl ./${1}.mpileup.out ./mRNA.bed > ./${1}.TranscriptionError
perl ATCG.pl  ${1}.mpileup.out 100 mRNA.bed > ${1}.ATGC20
