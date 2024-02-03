# GeneFuse vRust
This program is a Rust porting version based on [**GeneFuse 0.8.0**](https://github.com/OpenGene/GeneFuse) (original language: C++)


## Build from source
1. If you don't have Rust, [**install it**](https://www.rust-lang.org/learn/get-started).  
2. Build it with the below command on the root dir of this repo:
```
cargo build --release
```
The command will make a binary in {repo_root}/target/release/genefuse

3. Run genefuse
```
target/release/genefuse
```

## Notes: Some modifications in this version.
- Added read name as a criteria to sort `Match` objects. (The original codes can yield different unique read count per run.)
- It can accept **a file having a list of fusion csvs**. If you give it a file like that instead of a single csv, it outputs report files per fusion csv file.
- Parallelized part of `matcher::makeIndex()` method to increase performance.
- Some multi-threading codes were modified to be used in Rust.

## Performance test
```
genefuse \
    -r hg38.fa \
    -f genes/druggable.hg38.csv \
    -1 genefuse.R1.fq \
    -2 genefuse.R2.fq \
    -h report.html \
    -j report.json \
    -t 4
```
In a tiny test which ran genefuse 10 times with the above options,  
This version's speed was about **4.5x faster** on my machine. \[Rust:159s, Original:709s\] (WSL2 Ubuntu, Ryzen 5800x).


#### For more information about this program, please see [**the original repository**](https://github.com/OpenGene/GeneFuse).




