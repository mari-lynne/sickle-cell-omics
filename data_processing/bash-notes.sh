#!/bin/bash

srv=/fh/scratch/delete90/kooperberg_c/topmed_freeze#/minDP
wd=~/Documents/whi_sca/genomics

/run/user/1000/gvfs/sftp:host=rhino,user=mjohnso5/home/mjohnso5

sftp://mjohnso5@rhino/home/sct

plink2 --bcf ${srv}TOPMED38.vcf.gz --make-pfil--out ${wd}/TOP

#copied from sftp://mjohnso5@rhino/fh/scratch/delete90/koop_c/imputed10000g
plink2 --vcf WHI_share_Affy6.0-2015-03-05.chr11.vcf.gz --make-pgen--out whi_share_chr11
# then look in psam file for ids
# make id keep list

plink2 --pfile whi_share_chr11 --keep keep_list.txt --make-pgen --out whi_share_rna

# filter for sca variant

# Chromosome	11 Position	5227002

--snp chr:11