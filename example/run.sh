#!/bin/bash
perl ../muprimer -i ~/references/B73Ref4/genome/B73Ref4.fa \
	-g /homes/liu3zhen/references/B73Ref4/genemodel2/Zea_mays.B73_RefGen_v4.46.gtf \
	-m 1i_mu.txt -p out
