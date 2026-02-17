# Hyperemesis Gravidarum in MoBa

This repository contains code for genetic analyses of hyperemesis gravidarum in the [Norwegian Mother, Father and Child Cohort Study (MoBa)](https://www.fhi.no/en/ch/studies/moba). Genome-wide analyses were conducted on the [PsychGen release](https://doi.org/10.1101/2022.06.23.496289) of the genotypes. These analyses were in part included in the GWAS meta-analysis entitled _Multi-ancestry GWAS of severe pregnancy nausea and vomiting identifies risk loci associated with appetite, insulin signaling, and brain plasticity_. GWAS summary statistics are available at [DOI: 10.5281/zenodo.18274563](https://doi.org/10.5281/zenodo.18274563).

### Genome-wide association studies
The association of different phenotypes as conducted against the genomes of the children, mothers, and fathers in MoBa. GWAS summary statistics are available at [DOI: 10.5281/zenodo.18674365](https://doi.org/10.5281/zenodo.18674365).

- [Nausea vomiting during pregnancy](docs/23-07-21/nausea_vomiting.md): GWAS of participants of pregnancies where the mother reported suffering from nausea or vomiting versus mothers who did not report nausea or vomiting as control.
- [Hyperemesis gravidarum](docs/23-07-21/hyperemesis_gravidarum_vs_all.md): GWAS of participants of pregnancies where the mother reported being hospitalized due to prolonged nausea and vomiting versus all other pregnancies.
- [Hyperemesis gravidarum vs. no nausea vomiting](docs/23-07-21/hyperemesis_gravidarum_vs_no_nausea_vomiting.md): GWAS of participants of pregnancies where the mother reported being hospitalized due to prolonged nausea and vomiting versus mothers who did not report nausea or vomiting as control.
- [Nausea vomiting during pregnancy](docs/23-07-21/nausea_vomiting_strength.md): GWAS of participants of pregnancies encoding {0: mothers who did not report nausea or vomiting, 1: mothers reported nausea or vomiting, 2: mothers reported being hospitalized}.
- [Nausea before week 4](docs/23-07-21/nausea_before_4w.md): GWAS of participants of pregnancies where the mother reported nausea before week 4.
- [Vomiting before week 4](docs/23-07-21/vomiting_before_4w.md): GWAS of participants of pregnancies where the mother reported vomiting before week 4.
- [Nausea vomiting before week 4](docs/23-07-21/nausea_vomiting_before_4w.md): GWAS of participants of pregnancies where the mother reported nausea vomiting before week 4.
- [Nausea vomiting before week 8](docs/23-07-21/nausea_vomiting_before_8w.md): GWAS of participants of pregnancies where the mother reported nausea vomiting before week 8.
- [Nausea vomiting week 5 to 8](docs/23-07-21/nausea_vomiting_5w_8w.md): GWAS of participants of pregnancies where the mother reported nausea vomiting stratified by week.
- [Nausea vomiting week 9 to 12](docs/23-07-21/nausea_vomiting_9w_12w.md): GWAS of participants of pregnancies where the mother reported nausea vomiting stratified by week.
- [Nausea vomiting week 13 to 15](docs/23-07-21/nausea_vomiting_13w_15w.md): GWAS of participants of pregnancies where the mother reported nausea vomiting stratified by week.
- [Nausea vomiting week 13 to 16](docs/23-07-21/long_term_nausea_vomiting_13w_16w.md): GWAS of participants of pregnancies where the mother reported nausea vomiting stratified by week.
- [Nausea vomiting week 17 to 20](docs/23-07-21/long_term_nausea_vomiting_17w_20w.md): GWAS of participants of pregnancies where the mother reported nausea vomiting stratified by week.
- [Nausea vomiting week 21 to 24](docs/23-07-21/long_term_nausea_vomiting_21w_24w.md): GWAS of participants of pregnancies where the mother reported nausea vomiting stratified by week.
- [Nausea vomiting week 25 to 29](docs/23-07-21/long_term_nausea_vomiting_25w_28w.md): GWAS of participants of pregnancies where the mother reported nausea vomiting stratified by week.
- [Long-term nausea vomiting after week 12](docs/23-07-21/long_term_nausea_vomiting_after_13w.md): GWAS of participants of pregnancies where the mother reported suffering from long-term nausea vomiting after week 13.
- [Long-term nausea vomiting after week 17](docs/23-07-21/long_term_nausea_vomiting_after_17w.md): GWAS of participants of pregnancies where the mother reported suffering from long-term nausea vomiting after week 17.
- [Long-term nausea vomiting after week 21](docs/23-07-21/long_term_nausea_vomiting_after_21w.md): GWAS of participants of pregnancies where the mother reported suffering from long-term nausea vomiting after week 21.
- [Long-term nausea vomiting after week 25](docs/23-07-21/long_term_nausea_vomiting_after_25w.md): GWAS of participants of pregnancies where the mother reported suffering from long-term nausea vomiting after week 25.
- [Long-term nausea vomiting after week 29](docs/23-07-21/long_term_nausea_vomiting_after_29w.md): GWAS of participants of pregnancies where the mother reported suffering from long-term nausea vomiting after week 29.
- [Vomiting start](docs/23-07-21/vomiting_week_from.md): GWAS of participants of pregnancies against the week where mothers reported beginning to suffer from vomiting.
- [Vomiting end](docs/23-07-21/vomiting_week_to.md): GWAS of participants of pregnancies against the week where mothers reported ending to suffer from vomiting.
- [Vomiting duration](docs/23-07-21/vomiting_duration.md): GWAS of participants of pregnancies against the duration in weeks where mothers reported suffering from vomiting.

### Code
- [Genome-wide association study](src): Genome-wide association studies using regenie and top hit selection using cojo was run using the code in this folder.
- [WLM](src/look_up/utils/merge_and_wlm.R): WLM for top hits from the meta-analysis and their proxies was conducted using this script.
- [Other phenotypes](other_phenotypes/other_phenotypes.qmd): Association with other phenotypes was conducted using this notebook, results are exported in this [folder](other_phenotypes/tables).

#### Errors, questions, and bug report
We welcome bug reports, suggestions of improvements, and contributions. Please do not hesitate to open an issue a pull request in the repository.

#### Code of Conduct
As part of our efforts toward delivering open and inclusive science, we follow the [Contributor Convenant](https://www.contributor-covenant.org/) [Code of Conduct for Open Source Projects](CODE_OF_CONDUCT.md).

#### License
Unless otherwise specified in specific files, this work, the associated source code, and results are released under a [CC BY 4.0 License](https://creativecommons.org/licenses/by/4.0/).

![CC BY 4.0 License Logo](https://i.creativecommons.org/l/by/4.0/88x31.png)


