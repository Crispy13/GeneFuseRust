use std::sync::{Mutex, RwLock};

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use crossbeam::queue::ArrayQueue;

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

criterion_group!(benches, rwlock, mutex);
criterion_main!(benches);
