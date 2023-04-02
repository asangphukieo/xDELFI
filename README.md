# Nextflow Pipeline to generate DELFI and xDELFI feature
## Installation
1. Install docker (follows https://github.com/asangphukieo/CMUTEAM_NF_pipeline)
2. Install Nextflow (e.g. by anaconda: conda install -c bioconda nextflow)
3. Clone this repository https://github.com/asangphukieo/xDELFI
4. Corect paths in test.config
5. run ```nextflow run run_delfi.nf -with-report -with-dag dag.png -with-docker asangphukieo/xdelfi:v1.0.3 -c test.config -profile docker_test```
