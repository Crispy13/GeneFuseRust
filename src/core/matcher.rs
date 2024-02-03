use core::fmt;
use std::{
    collections::{hash_map::RandomState, BTreeMap, HashMap},
    mem,
    sync::{
        atomic::{AtomicI32, Ordering},
        Mutex,
    },
};

use genefuse::aux::int_hasher::{CPPTrivialHasherBuilder, FxHasherBuilder};
use rayon::iter::{IntoParallelRefIterator, ParallelBridge, ParallelIterator};

use super::{common::GenePos, fasta_reader::FastaReader, sequence::Sequence};

pub(crate) type GFHasherBuilder = FxHasherBuilder;

pub(crate) struct MatchResult {
    start_gp: GenePos,
    sequence: Sequence,
    reversed: bool,
    mismatches: Vec<i32>,
}

impl MatchResult {
    fn print(&self) {
        todo!()
    }
}

// we use 512M memory
const BLOOM_FILTER_LENGTH: usize = 1 << 29;
const KMER: i32 = 16;

pub(crate) struct Matcher<'f> {
    pub(crate) m_kmer_positions: HashMap<u32, Vec<GenePos>, GFHasherBuilder>,
    pub(crate) m_contig_names: Vec<String>,

    m_reference: Option<&'f mut FastaReader>,
    m_unique_pos: i32,
    m_dupe_pos: i32,
    m_bloom_filter_array: Box<[u8]>,
}

impl<'f> Matcher<'f> {
    pub(crate) fn from_ref_and_seqs(
        fasta_ref: Option<&'f mut FastaReader>,
        seqs: &[Sequence],
    ) -> Self {
        let mut matcher = Self {
            m_kmer_positions: HashMap::with_hasher(GFHasherBuilder::new()),
            m_contig_names: Vec::new(),
            m_reference: fasta_ref,
            m_unique_pos: 0,
            m_dupe_pos: 0,
            m_bloom_filter_array: vec![0; BLOOM_FILTER_LENGTH].into_boxed_slice(),
        };

        matcher.init_bloom_filter(seqs);
        matcher.make_index();

        matcher
    }

    fn init_bloom_filter(&mut self, seqs: &[Sequence]) -> () {
        let m_bloom_filter_array: [u8; BLOOM_FILTER_LENGTH] = [0; BLOOM_FILTER_LENGTH];

        for seq in seqs.iter() {
            self.init_bloom_filter_with_seq(seq);
            self.init_bloom_filter_with_seq(&seq.reverse_complement());
        }
    }

    fn init_bloom_filter_with_seq(&mut self, seq: &Sequence) {
        let s = seq.m_str.as_str();
        let mut valid = false;

        for i in (0..(s.len() - KMER as usize + 1)) {
            let kmer = make_kmer(s, i, &mut valid);
            if !valid {
                continue;
            }
            // set bloom filter
            *self
                .m_bloom_filter_array
                .get_mut(kmer.wrapping_shr(3) as usize)
                .unwrap() |= (1_i32.wrapping_shl((kmer & 0x07) as u32)) as u8;
        }
    }

    // fn make_index(&mut self) {
    //     if self.m_reference.is_none() {
    //         return;
    //     }

    //     let contig_ref = mem::take(&mut self.m_reference.as_mut().unwrap().m_all_contigs);

    //     let mut ctg = 0;

    //     log::debug!("indexing contig per ctg_name...");
    //     let mut seq_cv= Vec::new();
    //     for (ctg_name, s) in contig_ref.iter().map(|e| (e.0.as_str(), e.1.as_str())) {
    //         seq_cv.extend(s.as_bytes().iter().map(|b| b.to_ascii_uppercase()));
    //         // let s = s.to_uppercase();
    //         self.m_contig_names.push(ctg_name.to_owned());

    //         //index forward
    //         self.index_contig_bytes(ctg, &seq_cv, 0);

    //         //index reverse complement
    //         ctg += 1;
    //         seq_cv.clear();
    //     }

    //     log::info!("matcher indexing done");

    //     self.m_reference.as_mut().unwrap().m_all_contigs = contig_ref;

    // }

    fn make_index(&mut self) {
        if self.m_reference.is_none() {
            return;
        }

        let contig_ref = mem::take(&mut self.m_reference.as_mut().unwrap().m_all_contigs);
        let m_contig_names = Mutex::new(mem::take(&mut self.m_contig_names));
        let m_kmer_positions = Mutex::new(mem::take(&mut self.m_kmer_positions));

        // let ctg = AtomicI32::new(0);

        log::debug!("indexing contig per ctg_name...");
        // let mut seq_cv= Vec::new();
        contig_ref.iter().enumerate().par_bridge().for_each(|(ctg, e)| {
            let (ctg_name, s) = (e.0.as_str(), e.1.as_str());
            let seq_cv = s
                .as_bytes()
                .iter()
                .copied()
                .map(|b| b.to_ascii_uppercase())
                .collect::<Vec<u8>>();
            // let s = s.to_uppercase();
            m_contig_names.lock().unwrap().push(ctg_name.to_owned());

            //index forward
            self.index_contig_bytes(ctg as i32, &seq_cv, 0, &m_kmer_positions);

            //index reverse complement
            // ctg.fetch_add(1, Ordering::Relaxed);
            // seq_cv.clear();
        });

        log::info!("matcher indexing done");

        self.m_reference.as_mut().unwrap().m_all_contigs = contig_ref;
        self.m_kmer_positions = m_kmer_positions.into_inner().unwrap();
        self.m_contig_names = m_contig_names.into_inner().unwrap();
    }

    fn index_contig_cv(&mut self, ctg: i32, seq: &[char], start: i32) {
        let mut kmer = 0_u32;
        let mut valid = false;

        for (base, i) in seq
            .get(..((seq.len() as i32 - KMER) as usize))
            .unwrap()
            .iter()
            .zip(0..)
        {
            if valid {
                // let base = seq.get((i + KMER - 1) as usize).unwrap();
                let num = base2num(*base);

                if num < 0 {
                    valid = false;
                    continue;
                } else {
                    kmer = kmer.wrapping_shl(2) | num as u32;
                }
            } else {
                kmer = make_kmer_cv(seq, i as usize, &mut valid);
                if !valid {
                    continue;
                }
            }
            // check bloom filter
            if self
                .m_bloom_filter_array
                .get(kmer.wrapping_shr(3) as usize)
                .unwrap_or_else(|| {
                    panic!(
                        "kmer={}, idx={} bfa_len={}",
                        kmer,
                        kmer.wrapping_shr(3) as usize,
                        self.m_bloom_filter_array.len()
                    )
                })
                & (1_i32.wrapping_shl((kmer & 0x07) as u32)) as u8
                == 0
            {
                continue;
            }

            let site = GenePos {
                contig: ctg as i16,
                position: i + start,
            };

            self.m_kmer_positions
                .entry(kmer as u32)
                .or_insert_with(|| Vec::new())
                .push(site);
        }
    }

    fn index_contig_bytes(
        &self,
        ctg: i32,
        seq: &[u8],
        start: i32,
        kmer_positions: &Mutex<HashMap<u32, Vec<GenePos>, GFHasherBuilder>>,
    ) {
        let mut kmer = 0_u32;
        let mut valid = false;

        for (base, i) in seq
            .get(..((seq.len() as i32 - KMER) as usize))
            .unwrap()
            .iter()
            .zip(0..)
        {
            if valid {
                // let base = seq.get((i + KMER - 1) as usize).unwrap();
                let num = base2num_bytes(base);

                if num < 0 {
                    valid = false;
                    continue;
                } else {
                    kmer = kmer.wrapping_shl(2) | num as u32;
                }
            } else {
                kmer = make_kmer_bytes(seq, i as usize, &mut valid);
                if !valid {
                    continue;
                }
            }
            // check bloom filter
            if self
                .m_bloom_filter_array
                .get(kmer.wrapping_shr(3) as usize)
                .unwrap_or_else(|| {
                    panic!(
                        "kmer={}, idx={} bfa_len={}",
                        kmer,
                        kmer.wrapping_shr(3) as usize,
                        self.m_bloom_filter_array.len()
                    )
                })
                & (1_i32.wrapping_shl((kmer & 0x07) as u32)) as u8
                == 0
            {
                continue;
            }

            let site = GenePos {
                contig: ctg as i16,
                position: i + start,
            };

            kmer_positions
                .lock()
                .unwrap()
                .entry(kmer as u32)
                .or_insert_with(|| Vec::new())
                .push(site);
        }
    }

    // fn index_contig_bytes(&mut self, ctg: i32, seq: &[u8], start: i32) {
    //     let mut kmer = 0_u32;
    //     let mut valid = false;

    //     for (base, i) in seq.get(..((seq.len() as i32 - KMER) as usize)).unwrap().iter().zip(0..) {
    //         if valid {
    //             // let base = seq.get((i + KMER - 1) as usize).unwrap();
    //             let num = base2num_bytes(base);

    //             if num < 0 {
    //                 valid = false;
    //                 continue;
    //             } else {
    //                 kmer = kmer.wrapping_shl(2) | num as u32;
    //             }
    //         } else {
    //             kmer = make_kmer_bytes(seq, i as usize, &mut valid);
    //             if !valid {
    //                 continue;
    //             }
    //         }
    //         // check bloom filter
    //         if self
    //             .m_bloom_filter_array
    //             .get(kmer.wrapping_shr(3) as usize)
    //             .unwrap_or_else(||
    //                 panic!("kmer={}, idx={} bfa_len={}", kmer, kmer.wrapping_shr(3) as usize, self.m_bloom_filter_array.len())
    //             )
    //             & (1_i32.wrapping_shl((kmer & 0x07) as u32)) as u8
    //             == 0
    //         {
    //             continue;
    //         }

    //         let site = GenePos {
    //             contig: ctg as i16,
    //             position: i + start,
    //         };

    //         self.m_kmer_positions
    //             .entry(kmer as u32)
    //             .or_insert_with(|| Vec::new())
    //             .push(site);
    //     }
    // }

    fn index_contig(&mut self, ctg: i32, seq: &str, start: i32) {
        let mut kmer = 0_u32;
        let mut valid = false;

        for i in (0..(seq.len() - KMER as usize)).map(|e| e as i32) {
            if valid {
                let base = seq.chars().nth((i + KMER - 1) as usize).unwrap();
                let num = base2num(base);

                if num < 0 {
                    valid = false;
                    continue;
                } else {
                    kmer = kmer.wrapping_shl(2) | num as u32;
                }
            } else {
                kmer = make_kmer(seq, i as usize, &mut valid);
                if !valid {
                    continue;
                }
            }
            // check bloom filter
            if self
                .m_bloom_filter_array
                .get(kmer.wrapping_shr(3) as usize)
                .unwrap_or_else(|| {
                    panic!(
                        "kmer={}, idx={} bfa_len={}",
                        kmer,
                        kmer.wrapping_shr(3) as usize,
                        self.m_bloom_filter_array.len()
                    )
                })
                & (1_i32.wrapping_shl((kmer & 0x07) as u32)) as u8
                == 0
            {
                continue;
            }

            let site = GenePos {
                contig: ctg as i16,
                position: i + start,
            };

            self.m_kmer_positions
                .entry(kmer)
                .or_insert_with(|| Vec::new())
                .push(site);
        }
    }

    fn map_to_index(&mut self, sequence: &Sequence) -> Option<MatchResult> {
        let mut kmer_stat: HashMap<i64, i32, GFHasherBuilder> =
            HashMap::with_hasher(GFHasherBuilder::new());

        kmer_stat.insert(0, 0);

        let seq = sequence.m_str.as_str();

        const step: i32 = 1;
        const skip_threshold: i32 = 50;

        let seq_len = seq.len();

        let mut all_kmer = vec![0; seq_len];
        let mut kmer_valid = vec![false; seq_len];
        let mut skipped = vec![false; seq_len];

        // first pass, we only want to find if this seq can be partially aligned to the target
        let mut valid = false;
        let seq_bytes = seq.as_bytes();
        for i in (0..((seq_len as i32 - KMER + 1) as usize)) {
            let kmer = make_kmer_bytes(seq_bytes, i as usize, &mut valid);

            *kmer_valid.get_mut(i).unwrap() = valid;
            if !valid {
                continue;
            }

            *all_kmer.get_mut(i).unwrap() = kmer;
            // no match

            let kmer_pos_opt = self.m_kmer_positions.get(&kmer);

            if kmer_pos_opt.is_none() {
                *kmer_stat.get_mut(&0).unwrap() += 1;
                continue;
            }

            if kmer_pos_opt.unwrap().len() as i32 > skip_threshold {
                *skipped.get_mut(i).unwrap() = true;
                continue;
            }

            let kmer_pos = self.m_kmer_positions.get(&(kmer)).unwrap();
            for (i, gp) in kmer_pos.iter().enumerate() {
                let gp_i64 = gp_to_i64(&shift(gp, i as i32));
                kmer_stat
                    .entry(gp_i64)
                    .and_modify(|v| *v += 1)
                    .or_insert_with(|| 1);

                // if let Some(s) = kmer_stat.get_mut(&gp_i64) {
                //     *s += 1;
                // } else {
                //     kmer_stat.insert(gp_i64, 1);
                // }
            }
        }

        // get top N
        const TOP: usize = 5;

        let mut topgp = [0_i64; TOP];
        let mut topcount = [0_i32; TOP];

        for (gp, count) in kmer_stat.iter().map(|(a, b)| (*a, *b)) {
            // no need to update the top N
            if gp == 0 || count <= *topcount.get(TOP - 1).unwrap() {
                continue;
            }

            // update the last one first
            *topgp.get_mut(TOP - 1).unwrap() = gp;
            *topcount.get_mut(TOP - 1).unwrap() = count;

            // compare with the rest ones
            for t in (0..=(TOP - 2)).rev() {
                if count > *topcount.get(t).unwrap() {
                    *topcount.get_mut(t + 1).unwrap() = *topcount.get(t).unwrap();
                    *topgp.get_mut(t + 1).unwrap() = *topgp.get(t).unwrap();
                    *topcount.get_mut(t).unwrap() = count;
                    *topgp.get_mut(t).unwrap() = gp;
                }
            }
        }

        for t in (0..TOP) {
            if *topcount.get(t).unwrap() == 0 {
                break;
            }
            let mut mask = vec![0_u8; seq_len];

            // make the mask
            let mut valid = false;
            for i in (0..(seq_len as i32 - KMER + 1)) {
                valid = *kmer_valid.get(i as usize).unwrap();
                let kmer = *all_kmer.get(i as usize).unwrap();

                if !valid || self.m_kmer_positions.contains_key(&kmer) {
                    continue;
                }

                if !skipped.get(i as usize).unwrap()
                    && self.m_kmer_positions.get(&(kmer)).unwrap().len() < 5
                {
                    for gp in self.m_kmer_positions.get(&(kmer)).unwrap().iter() {
                        let gp_i64 = gp_to_i64(&shift(gp, i));
                        if (gp_i64 - topgp.get(t).unwrap()).abs() <= 2 {
                            for m in (i..(seq_len as i32).min(i + KMER)) {
                                *mask.get_mut(m as usize).unwrap() = 1;
                            }
                        }
                    }
                } else {
                    // this is repetive kmer, better method using binary search
                    if self.is_consistent(*topgp.get(t).unwrap(), kmer, i, 2) {
                        for m in (i..(seq_len as i32).min(i + KMER)) {
                            *mask.get_mut(m as usize).unwrap() = 1;
                        }
                    }
                }
            }
            let mut mismatches = Vec::<i32>::new();
            for i in (0..seq_len) {
                if *mask.get(i).unwrap() == 0 {
                    mismatches.push(i as i32);
                }
            }
            if mismatches.len() < 10 {
                let mr = MatchResult {
                    sequence: sequence.clone(),
                    start_gp: i64_to_gp(*topgp.get(t).unwrap()),
                    mismatches: mismatches,
                    reversed: false,
                };

                return Some(mr);
            }
        }

        None
    }

    // fn map_to_index(&mut self, sequence: &Sequence) -> Option<MatchResult> {
    //     let mut kmer_stat: HashMap<i64, i32> = HashMap::new();
    //     kmer_stat.insert(0, 0);

    //     let seq = sequence.m_str.as_str();

    //     const step: i32 = 1;
    //     const skip_threshold: i32 = 50;

    //     let seq_len = seq.len();

    //     let mut all_kmer = vec![0; seq_len];
    //     let mut kmer_valid = vec![false; seq_len];
    //     let mut skipped = vec![false; seq_len];

    //     // first pass, we only want to find if this seq can be partially aligned to the target
    //     let mut valid = false;
    //     for i in (0..((seq_len as i32 - KMER + 1) as usize)) {
    //         let kmer = make_kmer(seq, i as usize, &mut valid);

    //         *kmer_valid.get_mut(i).unwrap() = valid;
    //         if !valid {
    //             continue;
    //         }

    //         *all_kmer.get_mut(i).unwrap() = kmer;
    //         // no match
    //         if self.m_kmer_positions.get(&(kmer as u32)).is_none() {
    //             *kmer_stat.get_mut(&0).unwrap() += 1;
    //             continue;
    //         }

    //         if self.m_kmer_positions.len() as i32 > skip_threshold {
    //             *skipped.get_mut(i).unwrap() = true;
    //             continue;
    //         }

    //         let kmer_pos = self.m_kmer_positions.get(&(kmer as u32)).unwrap();
    //         for (i, gp) in kmer_pos.iter().enumerate() {
    //             let gp_i64 = gp_to_i64(&shift(gp, i as i32));
    //             kmer_stat
    //                 .entry(gp_i64)
    //                 .and_modify(|v| *v += 1)
    //                 .or_insert_with(|| 1);
    //         }
    //     }

    //     // get top N
    //     const TOP: usize = 5;

    //     let mut topgp = [0_i64; TOP];
    //     let mut topcount = [0_i32; TOP];

    //     for (gp, count) in kmer_stat.iter().map(|(a, b)| (*a, *b)) {
    //         // no need to update the top N
    //         if gp == 0 || count <= *topcount.get(TOP - 1).unwrap() {
    //             continue;
    //         }

    //         // update the last one first
    //         *topgp.get_mut(TOP - 1).unwrap() = gp;
    //         *topcount.get_mut(TOP - 1).unwrap() = count;

    //         // compare with the rest ones
    //         for t in (0..=(TOP - 2)).rev() {
    //             if count > *topcount.get(t).unwrap() {
    //                 *topcount.get_mut(t + 1).unwrap() = *topcount.get(t).unwrap();
    //                 *topgp.get_mut(t + 1).unwrap() = *topgp.get(t).unwrap();
    //                 *topcount.get_mut(t).unwrap() = count;
    //                 *topgp.get_mut(t).unwrap() = gp;
    //             }
    //         }
    //     }

    //     for t in (0..TOP) {
    //         if *topcount.get(t).unwrap() == 0 {
    //             break;
    //         }
    //         let mut mask = vec![0_u8; seq_len];

    //         // make the mask
    //         let mut valid = false;
    //         for i in (0..(seq_len as i32 - KMER + 1)) {
    //             valid = *kmer_valid.get(i as usize).unwrap();
    //             let kmer = *all_kmer.get(i as usize).unwrap();

    //             if !valid || self.m_kmer_positions.get(&(kmer as u32)).is_none() {
    //                 continue;
    //             }

    //             if !skipped.get(i as usize).unwrap()
    //                 && self.m_kmer_positions.get(&(kmer as u32)).unwrap().len() < 5
    //             {
    //                 for gp in self.m_kmer_positions.get(&(kmer as u32)).unwrap().iter() {
    //                     let gp_i64 = gp_to_i64(&shift(gp, i));
    //                     if (gp_i64 - topgp.get(t).unwrap()).abs() <= 2 {
    //                         for m in (i..(seq_len as i32).min(i + KMER)) {
    //                             *mask.get_mut(m as usize).unwrap() = 1;
    //                         }
    //                     }
    //                 }
    //             } else {
    //                 // this is repetive kmer, better method using binary search
    //                 if self.is_consistent(*topgp.get(t).unwrap(), kmer, i, 2) {
    //                     for m in (i..(seq_len as i32).min(i + KMER)) {
    //                         *mask.get_mut(m as usize).unwrap() = 1;
    //                     }
    //                 }
    //             }
    //         }
    //         let mut mismatches = Vec::<i32>::new();
    //         for i in (0..seq_len) {
    //             if *mask.get(i).unwrap() == 0 {
    //                 mismatches.push(i as i32);
    //             }
    //         }
    //         if mismatches.len() < 10 {
    //             let mr = MatchResult {
    //                 sequence: sequence.clone(),
    //                 start_gp: i64_to_gp(*topgp.get(t).unwrap()),
    //                 mismatches: mismatches,
    //                 reversed: false,
    //             };

    //             return Some(mr);
    //         }
    //     }

    //     None
    // }

    pub(crate) fn do_match(&mut self, sequence: &Sequence) -> Option<MatchResult> {
        let rcseq = sequence.reverse_complement();

        // log::info!("map_to_index() sequence...");
        let mut mco = self.map_to_index(sequence);

        if let Some(ref mut mc) = mco {
            mc.reversed = false;
        }

        // log::info!("map_to_index() rcseq...");
        let mut rcmco = self.map_to_index(&rcseq);
        if let Some(ref mut rcmc) = rcmco {
            rcmc.reversed = true;
        }

        if mco.is_none() {
            return rcmco;
        } else if rcmco.is_none() {
            return mco;
        } else {
            if mco.as_ref().unwrap().mismatches.len() <= rcmco.as_ref().unwrap().mismatches.len() {
                return mco;
            } else {
                return rcmco;
            }
        }
    }

    fn is_consistent(&self, thisgp: i64, kmer: u32, seqpos: i32, threshold: i32) -> bool {
        let gps = self.m_kmer_positions.get(&(kmer as u32)).unwrap();
        // align by seqpos
        let target = shift(&i64_to_gp(thisgp), -seqpos);
        let size = gps.len();
        let mut left = 0;
        let mut right = size - 1;

        while left <= right {
            let center = (left + right) / 2;
            let center_pos = gps.get(center).unwrap();

            if center_pos.contig < target.contig {
                // go right
                left = center + 1;
            } else if center_pos.contig > target.contig {
                // go left
                right = center - 1;
            } else {
                if (center_pos.position - target.position).abs() <= threshold {
                    return true;
                }

                if center_pos.position < target.position {
                    // go right
                    left = center + 1;
                } else if center_pos.position > target.position {
                    // go left
                    right = center - 1;
                }
            }
        }

        false
    }
}
#[inline]
fn base2num_bytes(c: &u8) -> i32 {
    match c {
        b'A' => {
            return 0;
        }
        b'T' => {
            return 1;
        }
        b'C' => {
            return 2;
        }
        b'G' => {
            return 3;
        }
        _ => {
            return -1;
        }
    }
}

#[inline]
fn base2num(c: char) -> i32 {
    match c {
        'A' => {
            return 0;
        }
        'T' => {
            return 1;
        }
        'C' => {
            return 2;
        }
        'G' => {
            return 3;
        }
        _ => {
            return -1;
        }
    }
}

fn make_kmer_cv(seq: &[char], pos: usize, valid: &mut bool) -> u32 {
    let mut kmer = 0_u32;
    for (i, base) in seq
        .get(pos..((pos as i32 + KMER) as usize))
        .unwrap()
        .iter()
        .enumerate()
    {
        match base {
            'A' => {
                kmer += 0;
                break;
            }
            'T' => {
                kmer += 1;
                break;
            }
            'C' => {
                kmer += 2;
                break;
            }
            'G' => {
                kmer += 3;
                break;
            }
            _ => {
                *valid = false;
                return 0;
            }
        }

        // not the tail
        if (i as i32) < (KMER - 1) {
            kmer = kmer.wrapping_shl(2);
        }
    }

    *valid = true;
    kmer
}

fn make_kmer_bytes(seq: &[u8], pos: usize, valid: &mut bool) -> u32 {
    let mut kmer = 0_u32;
    for (i, base) in seq
        .get(pos..((pos as i32 + KMER) as usize))
        .unwrap()
        .iter()
        .enumerate()
    {
        match base {
            b'A' => {
                kmer += 0;
                break;
            }
            b'T' => {
                kmer += 1;
                break;
            }
            b'C' => {
                kmer += 2;
                break;
            }
            b'G' => {
                kmer += 3;
                break;
            }
            _ => {
                *valid = false;
                return 0;
            }
        }

        // not the tail
        if (i as i32) < (KMER - 1) {
            kmer = kmer.wrapping_shl(2);
        }
    }

    *valid = true;
    kmer
}

fn make_kmer(seq: &str, pos: usize, valid: &mut bool) -> u32 {
    let mut kmer = 0_u32;
    for (i, base) in seq.chars().skip(pos).take(KMER as usize).enumerate() {
        match base {
            'A' => {
                kmer += 0;
                break;
            }
            'T' => {
                kmer += 1;
                break;
            }
            'C' => {
                kmer += 2;
                break;
            }
            'G' => {
                kmer += 3;
                break;
            }
            _ => {
                *valid = false;
                return 0;
            }
        }

        // not the tail
        if (i as i32) < (KMER - 1) {
            kmer = kmer.wrapping_shl(2);
        }
    }

    *valid = true;
    kmer
}

#[inline]
pub(crate) fn shift(gp: &GenePos, i: i32) -> GenePos {
    GenePos {
        contig: gp.contig,
        position: gp.position - i,
    }
}

#[inline]
pub(crate) fn gp_to_i64(gp: &GenePos) -> i64 {
    let ret = gp.contig as i64;

    (ret.wrapping_shl(32)) + gp.position as i64

    // gp.position
}

#[inline]
pub(crate) fn i64_to_gp(val: i64) -> GenePos {
    GenePos {
        contig: (val.wrapping_shr(32)) as i16,
        position: (val & 0x00000000FFFFFFFF) as i32,
    }
}
