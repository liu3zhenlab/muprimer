#!/bin/bash
ref=/homes/liu3zhen/references/B73Ref3/genome/Zea_mays.AGPv3.23.dna.genome.fa
gtf=/homes/liu3zhen/references/B73Ref3/genemodels/Zea_mays.AGPv3.27.gtf
mu=1i_mu.v3.txt
out=v3out
perl ../muprimer -i $ref -g $gtf -m $mu -p $out
