library(ensemblVEP)
# vep -i homo_sapiens_GRCh37.vcf --cache --port 3337 --buffer-size 8000 -o result
# less x*.vcf | paralle "vep --cache --port 3337 --force_overwrite -i {} -o vep_results/{.}.txt"