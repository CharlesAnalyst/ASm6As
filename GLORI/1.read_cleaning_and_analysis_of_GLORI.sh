trim_galore -q 20 --stringency 1 -e 0.3 --length 35
seqkit rmdup -s
fastx_trimmer -f 11

python run_GLORI.py -i GLORI-tools-main -q filter.fq -T 10 -f hg38.AG_conversion.fa -f2 hg38.fa -rvs hg38.rvsCom.fa -Tf GCF_000001405.40_GRCh38.p14_rna2.fna.AG_conversion.fa -a GCF_000001405.40_GRCh38.p14_genomic.gtf_change2Ens.tbl2 -b GCF_000001405.40_GRCh38.p14_genomic.gtf_change2Ens.tbl2.noredundance.base -pre sample_name -o out_dir --combine --rvs_fac
