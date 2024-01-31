use std::{collections::{BTreeMap, HashMap}, hash::BuildHasher, sync::{Mutex, RwLock}};

use criterion::{black_box, criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use crossbeam::queue::ArrayQueue;
use genefuse::aux::int_hasher::{CPPTrivialHasherBuilder, FxHasherBuilder};

use rustc_hash::FxHashMap;


#[derive(Default)]
pub struct CPPTrivialHasherBuilder2 {}

impl CPPTrivialHasherBuilder2 {
    pub fn new() -> Self {
        Self {}
    }
}

impl BuildHasher for CPPTrivialHasherBuilder2 {
    type Hasher = CPPTrivialHasher2;

    fn build_hasher(&self) -> Self::Hasher {
        CPPTrivialHasher2 { state: 0 }
    }
}

#[derive(Default)]
pub struct CPPTrivialHasher2 {
    state: u64,
}

impl std::hash::Hasher for CPPTrivialHasher2 {
    fn finish(&self) -> u64 {
        self.state
    }

    fn write(&mut self, bytes: &[u8]) {
        unreachable!()
    }

    fn write_i64(&mut self, i: i64) {
        // println!("input for hashing: {i}");
        self.state = u64::from_ne_bytes(i.to_ne_bytes());
        // println!("state= {}", self.state);
    }

    fn write_u32(&mut self, i: u32) {
        self.state = i as u64;  
    }
}


fn _rwlock(i: &(ArrayQueue<u8>, RwLock<u8>)) {
    let (aq, rl) = i;

    let lock = rl.write().unwrap();

    aq.push(0).unwrap();

    aq.pop().unwrap();

    drop(lock);
}

pub fn rwlock(c: &mut Criterion) {
    let aq = ArrayQueue::<u8>::new(1000);
    let rl = RwLock::<u8>::new(0);

    let i = (aq, rl);
    c.bench_with_input(BenchmarkId::new("rwlock", 1), &i, |b, i| {
        b.iter(|| _rwlock(i));
    });

    // c.bench_function("rwlock", |b| b.iter(|| fibonacci(black_box(20))));
}

fn _mutex(i: &(ArrayQueue<u8>, Mutex<u8>)) {
    let (aq, rl) = i;

    let lock = rl.lock().unwrap();

    aq.push(0).unwrap();

    aq.pop().unwrap();

    drop(lock);
}

pub fn mutex(c: &mut Criterion) {
    let aq = ArrayQueue::<u8>::new(1000);
    let rl = Mutex::<u8>::new(0);

    let i = (aq, rl);
    c.bench_with_input(BenchmarkId::new("mutex", 2), &i, |b, i| {
        b.iter(|| _mutex(i));
    });

    // c.bench_function("rwlock", |b| b.iter(|| fibonacci(black_box(20))));
}



const REPEAT:usize = 10000;
const START:usize = REPEAT / 2;
const END:usize = REPEAT * 3 / 4;
const TAKE:usize = END-START;

pub fn str_rev(c: &mut Criterion) {
    let a = "ATISJDFIWEOSLDKFJLKSDJFASOIDJAS".repeat(REPEAT);

    c.bench_function("str_rev", |b| {
        b.iter_batched(|| a.clone(), |s| {
            let mut s = s.clone().into_bytes();
            s.reverse();
            String::from_utf8(s).unwrap();
        },
    BatchSize::LargeInput);
    });

    // c.bench_function("rwlock", |b| b.iter(|| fibonacci(black_box(20))));
}

pub fn chars_rev(c: &mut Criterion) {
    let a = "ATISJDFIWEOSLDKFJLKSDJFASOIDJAS".repeat(REPEAT);

    c.bench_function("chars_rev", |b| {
        b.iter_batched(|| a.clone(), |s| {
            s.chars().rev().collect::<String>()
        },
    BatchSize::LargeInput);
    });

    // c.bench_function("rwlock", |b| b.iter(|| fibonacci(black_box(20))));
}


pub fn trivial_hash(c: &mut Criterion) {
    let a = (-25000..25000_i64).collect::<Vec<_>>();

    c.bench_function("trivial_hash", |b| {
        b.iter_batched(|| a.as_slice(), |s| {
            let mut map = HashMap::with_hasher(CPPTrivialHasherBuilder {});

            a.iter().for_each(|e| {map.insert(*e, 100);});
            
        },
    BatchSize::SmallInput);
    });

    // c.bench_function("rwlock", |b| b.iter(|| fibonacci(black_box(20))));
}

pub fn default_hash(c: &mut Criterion) {
    let a = (-25000..25000_i64).collect::<Vec<_>>();

    c.bench_function("default_hash", |b| {
        b.iter_batched(|| a.as_slice(), |s| {
            let mut map = HashMap::new();

            a.iter().for_each(|e| {map.insert(*e, 100);});
            
        },
    BatchSize::SmallInput);
    });

    // c.bench_function("rwlock", |b| b.iter(|| fibonacci(black_box(20))));
}

pub fn fx_hash(c: &mut Criterion) {
    let a = (-25000..25000_i64).collect::<Vec<_>>();

    c.bench_function("fx_hash", |b| {
        b.iter_batched(|| a.as_slice(), |s| {
            let mut map = HashMap::with_hasher(FxHasherBuilder::default());

            a.iter().for_each(|e| {map.insert(*e, 100);});
            
        },
    BatchSize::SmallInput);
    });

    // c.bench_function("rwlock", |b| b.iter(|| fibonacci(black_box(20))));
}


pub fn trivial_hash2(c: &mut Criterion) {
    let a = (-25000..25000_i64).collect::<Vec<_>>();

    c.bench_function("trivial_hash2", |b| {
        b.iter_batched(|| a.as_slice(), |s| {
            let mut map = HashMap::with_hasher(CPPTrivialHasherBuilder2 {});

            a.iter().for_each(|e| {map.insert(*e, 100);});
            
        },
    BatchSize::SmallInput);
    });

    // c.bench_function("rwlock", |b| b.iter(|| fibonacci(black_box(20))));
}

pub fn btree_hm(c: &mut Criterion) {
    let a = (-25000..25000_i64).collect::<Vec<_>>();

    c.bench_function("btree_hm", |b| {
        b.iter_batched(|| a.as_slice(), |s| {
            let mut map = BTreeMap::default();

            a.iter().for_each(|e| {map.insert(*e, 100);});
            
        },
    BatchSize::SmallInput);
    });

    // c.bench_function("rwlock", |b| b.iter(|| fibonacci(black_box(20))));
}

criterion_group!(benches, default_hash, trivial_hash, fx_hash, trivial_hash2, btree_hm);
criterion_main!(benches);
