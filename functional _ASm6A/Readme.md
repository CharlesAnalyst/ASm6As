# Input

`ASm6A_dir` contain files below :
The file (tab character delimitation) include 6 columns :
(produced by step2, combined significant and insignificant variants which overlap with m6A peak)

1. **contig** - The name of the chromosome.
2. **start** - The start position of the variant in the chromosome. The first base in a chromosome is numbered 0.
3. **position** -  The end position of the variant in the chromosome. The first base in a chromosome is numbered 0.
4. **term** - The information combined from refAllele, altAllele, allelicRatio and allelic direction of this variant.
5. **refRPKM_ratio** - m6A level in reference allele = ref_RPKM_IP/ref_RPKM_Input
6. **altRPKM_ratio** - m6A level in alternative allele = alt_RPKM_IP/alt_RPKM_Input
7. **mark** - The allelic direction of a variant. (When the variant is not statistically significant, this column will be blank.)
```
chr9	90321582	90321583	C>T;0.23;alt	3.5	11.53	alt
chr9    97849089	97849090	A>G;0.61;ref	2.32	1.49	ref
chr1	20977448	20977449	G>T;0.51;unsig	5.75	5.47	
chr1	55207537	55207538	G>A;0.47;unsig	6.16	6.84
```

`result_dir`  The directory which saved functional ASm6A file

`prop`  The threshold used in article to filter the ASm6A in same direction across samples

# Output

The output file  (tab character delimitation) include 3 columns :
1. **contig** - The name of the chromosome.
2. **start** - The start position of the variant in the chromosome. The first base in a chromosome is numbered 0.
3. **position** -  The end position of the variant in the chromosome. The first base in a chromosome is numbered 0.
```text
chr15	86288654	86288655
chr16	57116457	57116458
chr16	72110322	72110323
```