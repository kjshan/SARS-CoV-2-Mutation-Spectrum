#data from Table S3 from the paper below
#Garcia-Nieto, P.E., Morrison, A.J., and Fraser, H.B. (2019). The somatic mutation landscape of the human body. Genome Biol 20, 298.

#call somatic mutations from 36 tissues
#Somatic mutations were polarized according to the coding strand DNA based on the human genome annotation (Ensembl, GRCh37, version 84)
perl /Dell/Dell13/shankj/projects/Cov/SomaticMutations/Somatic.pl TableS3.tsv > Somatic.txt
