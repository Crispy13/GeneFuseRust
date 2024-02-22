# GeneFuse vRust
A Rust ported GeneFuse with improved performance.  

This program is a Rust porting version based on [**GeneFuse 0.8.0**](https://github.com/OpenGene/GeneFuse)


## How to use: Build from source
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
<br>

## Notes: Some modifications in this version
- Added read name as a criteria to sort `Match` objects. (The original codes can yield different unique read count per run.)
- It can accept **a file having a list of fusion csvs**. If you give it a file like that instead of a single csv, it outputs report files per fusion csv file.
- Parallelized part of `matcher::makeIndex()` method to increase performance.
- Some multi-threading codes were modified to be used in Rust.


## Performance test
* In a test, this version's running time was ~7.4x and it used 105% memory.
* [**Details**](./benchmark_res/bench_res.md)


## Multiple csv path file Example
```
./EGFR_KRAS.csv
./TP53_EGFR.csv
```


#### For more information about this program, please see [**the original repository**](https://github.com/OpenGene/GeneFuse).




