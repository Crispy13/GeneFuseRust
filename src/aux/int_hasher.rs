use std::{collections::HashMap, hash::BuildHasher};

use rustc_hash::FxHasher;

#[derive(Default)]
pub struct CPPTrivialHasherBuilder {}

impl CPPTrivialHasherBuilder {
    pub fn new() -> Self {
        Self {}
    }
}

impl BuildHasher for CPPTrivialHasherBuilder {
    type Hasher = CPPTrivialHasher;

    fn build_hasher(&self) -> Self::Hasher {
        CPPTrivialHasher { state: 0 }
    }
}

#[derive(Default)]
pub struct CPPTrivialHasher {
    state: u64,
}

impl std::hash::Hasher for CPPTrivialHasher {
    fn finish(&self) -> u64 {
        self.state
    }

    fn write(&mut self, bytes: &[u8]) {
        unreachable!()
    }

    fn write_i64(&mut self, i: i64) {
        // println!("input for hashing: {i}");
        self.state = i as u64;
        // println!("state= {}", self.state);
    }

    fn write_u32(&mut self, i: u32) {
        self.state = i as u64;  
    }
}


#[derive(Default)]
pub struct FxHasherBuilder {}

impl FxHasherBuilder {
    pub fn new() -> Self {
        Self {}
    }
}

impl BuildHasher for FxHasherBuilder {
    type Hasher = FxHasher;

    fn build_hasher(&self) -> Self::Hasher {
        FxHasher::default()
    }
}

#[cfg(test)]
mod test {
    use std::{collections::HashMap, ops::Neg};
    use super::*;

    #[test]
    fn htest() {
        let mut map = HashMap::<i64, i32, CPPTrivialHasherBuilder>::with_hasher(CPPTrivialHasherBuilder {});

        map.insert(1, 10);
        map.insert(-234234234, 10);

        // map.hasher().

        println!("{}", u64::from_le_bytes((234234234_i64).neg().to_le_bytes()));

    }
}


