##
# Plink commands used to extract the genotypes in an R-friendly format
# Plink2 run from conda
##

# Make pgen files with the ids and variants of interest
plink2 --bfile /mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/plink/19 --extract rare_variant --make-pgen --out /mnt/work/marc/hg/variants/rare_variant
plink2 --pfile /mnt/work/marc/hg/variants/rare_variant --export-allele rare_variant --export A --out /mnt/work/marc/hg/variants/rare_variant

