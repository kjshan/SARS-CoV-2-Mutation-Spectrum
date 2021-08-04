#ARC-seq
#call mRNA mutations and 20S RNA virus mutations

cd /Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/STAR/WT_Ctrl/
STAR  --runThreadN 30 --genomeDir Yeast_virus_STAR_index/ --readFilesIn  Ctrl_1.fq.gz Ctrl_2.fq.gz --outFileNamePrefix WT_Ctrl --outFilterType BySJout --readFilesCommand zcat

perl /Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/STAR/bin/Filter_Read.pl /Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/STAR/WT_Ctrl/WT_CtrlAligned.out.sam > WT_Ctrl.noIndel.Match.sam
samtools sort -@ 20 WT_Ctrl.noIndel.Match.sam -o WT_Ctrl.noIndel.Match.bam
samtools mpileup -Q 30 -d 0  --reference /Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ours/H2O2_Test/Clean/STAR/Yeast_virus_Genome.fasta  --output-QNAME  --output-BP  WT_Ctrl.noIndel.Match.bam -o  WT_Ctrl.noIndel.Match.mpileupQ30
perl /Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/STAR/bin/mpileup_out.pl WT_Ctrl.noIndel.Match.mpileupQ30 > WT_Ctrl.noIndel.Match.mpileup.outQ30
perl /Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/STAR/bin/mRNA_error.pl WT_Ctrl.noIndel.Match.mpileup.outQ30 /Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/Reference/mRNA.bed > WT_Ctrl.noIndel.Match.mRNA.errorQ30
perl /Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/STAR/bin/mRNA_base_count.pl WT_Ctrl.noIndel.Match.mpileup.outQ30 /Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/Reference/mRNA.bed > WT_Ctrl.100.Q30
