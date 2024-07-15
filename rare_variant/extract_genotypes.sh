##
# Plink commands used to extract the genotypes in an R-friendly format
# Plink2 run from conda
##

# Make pgen files with the ids and variants of interest
plink2 --bfile /mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/plink/19 --extract rare_variant_19 --make-pgen --out /mnt/work/marc/hg/variants/rare_variant_19
plink2 --pfile /mnt/work/marc/hg/variants/rare_variant_19 --export-allele rare_variant_19 --export A --out /mnt/work/marc/hg/variants/rare_variant_19
plink2 --bfile /mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/plink/10 --extract rare_variant_10 --make-pgen --out /mnt/work/marc/hg/variants/rare_variant_10
plink2 --pfile /mnt/work/marc/hg/variants/rare_variant_10 --export-allele rare_variant_10 --export A --out /mnt/work/marc/hg/variants/rare_variant_10

