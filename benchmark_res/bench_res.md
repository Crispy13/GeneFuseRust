### Bench 1: Paired End Scanning (vs [Original](https://github.com/OpenGene/GeneFuse))
|TestNum|Running Time(C)|Mem Usage(C)|Running Time(E)|Mem Usage(E)|Speed ratio(E to C)|Mem Usage Ratio(E to C)|
|---|---|---|---|---|---|---|
|#1|193.65s|8.1G|28.04s|8.13G|691%|100%|
|#2|193.75s|8.16G|28.02s|8.23G|691%|101%|
|#3|120.26s|5.02G|16.29s|5.62G|738%|112%|
|#4|120.34s|5.08G|18.64s|5.78G|646%|114%|
|#5|110.44s|6.06G|18.41s|5.72G|600%|94%|
|#6|79.57s|4.85G|13.11s|5.14G|607%|106%|

* C: Control, original version, E: Experimental, this version.
* Mean speed ratio = 662%
* Mean mem usage ratio = 105%



### Bench 1: Test options
```
tn="mc" #1
target/release/genefuse \
    -r genefuse_rust/hg19.fa \
    -f genefuse_rust/genes/cancer.hg19.csv \
    -1 genefuse_rust/genefuse.R1.fq \
    -2 genefuse_rust/genefuse.R2.fq \
    -h pe_report_${tn}.html \
    -j pe_report_${tn}.json

tn="mc38" #2
target/release/genefuse \
    -r genefuse_rust/hg38.fa \
    -f genefuse_rust/genes/cancer.hg38.csv \
    -1 genefuse_rust/genefuse.R1.fq \
    -2 genefuse_rust/genefuse.R2.fq \
    -h pe_report_${tn}.html \
    -j pe_report_${tn}.json

tn="md" #3
target/release/genefuse \
    -r genefuse_rust/hg19.fa \
    -f genefuse_rust/genes/druggable.hg19.csv \
    -1 genefuse_rust/genefuse.R1.fq \
    -2 genefuse_rust/genefuse.R2.fq \
    -h pe_report_${tn}.html \
    -j pe_report_${tn}.json

tn="md38" #4
target/release/genefuse \
    -r genefuse_rust/hg38.fa \
    -f genefuse_rust/genes/druggable.hg38.csv \
    -1 genefuse_rust/genefuse.R1.fq \
    -2 genefuse_rust/genefuse.R2.fq \
    -h pe_report_${tn}.html \
    -j pe_report_${tn}.json

tn="mtc" #5
target/release/genefuse \
    -r genefuse_rust/hg19.fa \
    -f genefuse_rust/testdata/cancer.csv \
    -1 genefuse_rust/testdata/R1.fq \
    -2 genefuse_rust/testdata/R2.fq \
    -h pe_report_${tn}.html \
    -j pe_report_${tn}.json

tn="mtf" #6
target/release/genefuse \
    -r genefuse_rust/hg19.fa \
    -f genefuse_rust/testdata/fusions.csv \
    -1 genefuse_rust/testdata/R1.fq \
    -2 genefuse_rust/testdata/R2.fq \
    -h pe_report_${tn}.html \
    -j pe_report_${tn}.json
```



### Bench 2: Multi csv input mode (vs [GeneFuse_Plus](https://github.com/tsy19900929/GeneFuse_Plus))
|TestNum|Running Time(C)|Mem Usage(C)|Running Time(E)|Mem Usage(E)|Speed ratio(E to C)|Mem Usage Ratio(E to C)|
|---|---|---|---|---|---|---|
|#1|6328.27s|9.59G|166.65s|16.65G|3797%|174%|

### Bench 2: Test Options
```
tn="hg38l"
target/release/genefuse \
    -r genefuse_rust/hg38.fa \
    -f hg38_fusion_csv_list.txt \
    -1 genefuse_rust/genefuse.R1.fq \
    -2 genefuse_rust/genefuse.R2.fq \
    -t 4 \
    -h mfmp_report_${tn}.html \
    -j mfmp_report_${tn}.json
```

### Test machine specs
```
System: Windows 11 WSL2 Ubuntu 22.04.4 LTS
CPU: AMD Ryzen 7 5800x 
MEM: 32GB
NVME SSD
```