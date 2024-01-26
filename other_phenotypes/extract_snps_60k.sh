##
# Plink commands used to extract the genotypes in an R-friendly format
# Plink2 run from conda
##

# Make pgen files with the ids and variants of interest
plink2 --bfile /mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc --extract /mnt/work/marc/github/moba/hyperemesis_gravidarum/src/look_up/resources/look_up_ids_23-07-21 --keep /mnt/archive/moba/pheno/v12/pheno_anthropometrics_24-01-16/id/mothers_id_unrelated_plink --make-pgen --out /mnt/work/marc/hg/variants/mother_genotypes

plink2 --bfile /mnt/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc --extract /mnt/work/marc/github/moba/hyperemesis_gravidarum/src/look_up/resources/look_up_ids_23-07-21 --keep /mnt/archive/moba/pheno/v12/pheno_anthropometrics_24-01-16/id/children_id_unrelated_plink --make-pgen --out /mnt/work/marc/hg/variants/children_genotypes

