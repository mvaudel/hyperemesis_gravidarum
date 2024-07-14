##
# Plink commands used to extract the genotypes in an R-friendly format
# Plink2 run from conda
##

# Make pgen files with the ids and variants of interest
plink2 --bfile /mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/plink/19 --extract rs372120002 --make-pgen --out /mnt/work/marc/hg/variants/rs372120002
plink2 --pfile /mnt/work/marc/hg/variants/rs372120002 --export-allele rs372120002 --export A --out /mnt/work/marc/hg/variants/rs372120002

