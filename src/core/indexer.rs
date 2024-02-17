// use super::fastq_reader::FastqReader;

use std::{
    array,
    collections::{hash_map, BTreeMap, HashMap},
    error,
    fmt::{self, Write},
    hash::Hash,
    io::Write as io_write,
    ops::Index, sync::Arc,
};

use crate::{
    aux::{
        global_settings::{global_settings, GlobalSettings},
        pbar::prepare_pbar,
    },
    core::{
        common::{DUPE_HIGH_LEVEL, DUPE_NORMAL_LEVEL},
        matcher::GFHasherBuilder,
        pescanner::DBT,
    },
    utils::{write_tsv_row, ObjectTsvWriter},
};

use super::{
    common::GenePos, fasta_reader::FastaReader, fastq_reader::FastqReader, fusion::Fusion, fusion_scan::Error, gene::Gene, read::SequenceRead, sequence::{reverse_complement, Sequence}
};

const MATCH_TOP: u8 = 3;
const MATCH_SECOND: u8 = 2;
const MATCH_NONE: u8 = 1;
const MATCH_UNKNOWN: u8 = 0;

const KMER: i32 = 16;
// 512M bloom filter
const BLOOM_FILTER_SIZE: usize = (1_i32.wrapping_shl(29)) as usize;
const BLOOM_FILTER_BITS: usize = BLOOM_FILTER_SIZE - 1;

#[derive(Debug)]
pub(crate) struct SeqMatch {
    pub(crate) seq_start: i32,
    pub(crate) seq_end: i32,
    pub(crate) start_gp: GenePos,
}

impl SeqMatch {
    fn new(seq_start: i32, seq_end: i32, start_gp: GenePos) -> Self {
        Self {
            seq_start,
            seq_end,
            start_gp,
        }
    }
}

impl fmt::Display for SeqMatch {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}-{}|{}:{}",
            self.seq_start, self.seq_end, self.start_gp.contig, self.start_gp.position
        )
    }
}

pub(crate) struct Indexer {
    m_ref_file: String,
    pub(crate) m_reference: Option<Arc<FastaReader>>,
    m_fusions: Vec<Fusion>,
    m_unique_pos: i32,
    m_dupe_pos: i32,

    pub(crate) m_kmer_pos: HashMap<i64, GenePos, GFHasherBuilder>,
    pub(crate) m_bloom_filter: Box<[u8]>, // `u8` in rust corresponds to unsigned char in C
    pub(crate) m_dupe_list: Vec<Vec<GenePos>>,
    pub(crate) m_fusion_seq: Vec<String>,
}

impl Indexer {
    pub(crate) fn new(ref_file: &str, fusions: Vec<Fusion>) -> Result<Self, Error> {
        let mut m_reference = FastaReader::new(&ref_file, false)?;

        log::debug!("Reading reference, {}", &ref_file);
        m_reference.read_all();

        Ok(Self {
            m_ref_file: ref_file.to_string(),
            m_reference: Some(Arc::new(m_reference)),
            m_fusions: fusions,
            m_unique_pos: 0,
            m_dupe_pos: 0,
            m_kmer_pos: HashMap::with_hasher(GFHasherBuilder::default()),
            m_bloom_filter: vec![0; BLOOM_FILTER_SIZE].into_boxed_slice(),
            m_dupe_list: Vec::new(),
            m_fusion_seq: Vec::new(),
        })
    }

    pub(crate) fn with_loaded_ref(m_reference: Arc<FastaReader>, fusions: Vec<Fusion>) -> Self {
        Self {
            m_ref_file: m_reference.m_fasta_file.clone(),
            m_reference: Some(m_reference),
            m_fusions: fusions,
            m_unique_pos: 0,
            m_dupe_pos: 0,
            m_kmer_pos: HashMap::with_hasher(GFHasherBuilder::default()),
            m_bloom_filter: vec![0; BLOOM_FILTER_SIZE].into_boxed_slice(),
            m_dupe_list: Vec::new(),
            m_fusion_seq: Vec::new(),
        }
    }

    pub(crate) fn get_ref(&self) -> Option<&FastaReader> {
        self.m_reference.as_ref().map(|e| e.as_ref())
    }

    // pub(crate) fn get_ref_mut(&mut self) -> Option<&mut FastaReader> {
    //     self.m_reference.as_mut()
    // }

    pub(crate) fn make_index(&mut self) {
        if self.m_reference.is_none() {
            return;
        }

        log::debug!("Index iter len={}", self.m_fusions.len());
        let pbar = prepare_pbar(self.m_fusions.len() as u64);
        pbar.set_message("making index...");

        //mutables : fusion seq
        for ctg in (0..self.m_fusions.len()).map(|e| {
            pbar.inc(1);
            e
        }) {
            let gene = &self.m_fusions.get(ctg).unwrap().m_gene;
            let mut chr = gene.m_chr.clone();

            // log::debug!("gene={:?}", gene);
            let seq_ref = &self.m_reference.as_ref().unwrap().m_all_contigs;
            if !seq_ref.contains_key(&chr) {
                if let Some(key) = check_map_has_key_then_get_back(seq_ref, format!("chr{}", chr)) {
                    chr = key;
                } else if let Some(key) =
                    check_map_has_key_then_get_back(seq_ref, chr.replace("chr", ""))
                {
                    chr = key;
                } else {
                    self.m_fusion_seq.push("".to_string());
                    continue;
                }
            }

            let s = seq_ref
                .get(&chr)
                .unwrap()
                .get((gene.m_start as usize)..(gene.m_end as usize))
                .unwrap()
                .to_uppercase();

            log::debug!("Indexing contig forward...",);
            //index forward
            self.index_contig(ctg, &s, 0);

            log::debug!("Indexing contig reverse...",);
            //index reverse complement
            let rev_comp_seq = reverse_complement(&s);
            self.index_contig(ctg, &rev_comp_seq, 1 - (s.len() as i32));

            self.m_fusion_seq.push(s);
        }

        pbar.finish_and_clear();

        self.fill_bloom_filter();
        log::info!("mapper indexing done.");
    }

    fn index_contig(&mut self, ctg: usize, seq: &str, start: i32) {
        let mut kmer = -1_i64;

        // let seq_cv = seq.chars().collect::<Vec<char>>();

        log::debug!(
            "Index contig={ctg}, start={start}, seq_char_len={}",
            seq.len()
        );
        for i in (0..(seq.len() as i32 - KMER)) {
            kmer = make_kmer_bytes(seq, i, kmer, 1);
            if kmer < 0 {
                continue;
            }

            let site = GenePos {
                contig: ctg as i16,
                position: i as i32 + start,
            };
            // site.contig = ctg as i16;
            // site.position = i as i32 + start;

            // this is a dupe
            if let Some(gp) = self.m_kmer_pos.get(&kmer).and_then(|gp| Some(gp.clone())) {
                // already marked as a dupe
                if gp.contig == DUPE_HIGH_LEVEL {
                    // skip it since it's already a high level dupe
                    continue;
                } else if gp.contig == DUPE_NORMAL_LEVEL {
                    if self.m_dupe_list.get(gp.position as usize).unwrap().len()
                        >= global_settings().skip_key_dup_threshold
                    {
                        self.m_kmer_pos.get_mut(&kmer).unwrap().contig = DUPE_HIGH_LEVEL;
                        *self.m_dupe_list.get_mut(gp.position as usize).unwrap() = Vec::new();
                    } else {
                        self.m_dupe_list
                            .get_mut(gp.position as usize)
                            .unwrap()
                            .push(site);
                    }
                } else {
                    // else make a new dupe entry
                    let mut gps = Vec::new();
                    gps.push(gp);
                    gps.push(site);
                    self.m_dupe_list.push(gps);

                    // and mark it as a dupe
                    {
                        let kmer_pos = self.m_kmer_pos.get_mut(&kmer).unwrap();
                        kmer_pos.contig = DUPE_NORMAL_LEVEL;
                        kmer_pos.position = (self.m_dupe_list.len() - 1) as i32;
                    }
                    self.m_unique_pos -= 1;
                    self.m_dupe_pos += 1;
                }
            } else {
                // log::debug!("kmer={}, m_kmer_pos={:#?}", kmer, self.m_kmer_pos);
                self.m_kmer_pos.insert(kmer, site);
                self.m_unique_pos += 1;
            }
        }
    }

    fn fill_bloom_filter(&mut self) -> () {
        for (kmer, gp) in self.m_kmer_pos.iter() {
            *self
                .m_bloom_filter
                .get_mut((kmer.wrapping_shr(3)) as usize)
                .unwrap() |= 0x1_i32.wrapping_shl((kmer & 0x07) as u32) as u8;
        }
    }

    pub(crate) fn map_read(&self, r: &SequenceRead) -> Vec<SeqMatch> {
        // log::debug!("m_kmer_pos={:#?}", self.m_kmer_pos);

        // let mut kmer_stat: HashMap<i64, i32, GFHasherBuilder> =
        //     HashMap::with_hasher(GFHasherBuilder::default());

        let mut kmer_stat: BTreeMap<i64, i32> = BTreeMap::new();
        kmer_stat.insert(0, 0);

        let seq = r.m_seq.m_str.as_str();
        const step: usize = 2_usize;
        let seq_cv = seq.chars().collect::<Vec<char>>();
        let seqlen = seq_cv.len() as i32;

        // let mut ocw = if r.m_name.contains(DBT) {
        //     let ocw = ObjectTsvWriter::default();

        //     Some(ocw)
        // } else {
        //     None
        // };

        // first pass, we only want to find if this seq can be partially aligned to the target
        let mut kmer = -1_i64;
        // log::debug!("map_read() -> seqlen={}", seqlen);
        for i in (0..(seqlen as i32 - KMER + 1)).step_by(step) {
            kmer = make_kmer_cv(&seq_cv, i, kmer, step as i32);
            // log::debug!("kmer={}, i={}", kmer, i);

            // ocw.as_mut().and_then(|mut f| {
            //     write_tsv_row!(&mut f, i,kmer, "1st-pass");
            //     Some(())
            // });

            if kmer < 0 {
                continue;
            }

            let pos = kmer.wrapping_shr(3);
            let bit = kmer & 0x07;
            if self.m_bloom_filter.get(pos as usize).unwrap()
                & (0x1_i64.wrapping_shl(bit as u32) as u8)
                == 0
            {
                *kmer_stat.get_mut(&0).unwrap() += 1;
                continue;
            }

            let gp = self.m_kmer_pos.get(&kmer).unwrap();
            // ocw.as_mut().and_then(|mut f| {
            //     write_tsv_row!(&mut f, i, gp.contig, "1st-pass");
            //     Some(())
            // });
            if gp.contig == DUPE_HIGH_LEVEL {
                // too much keys in this dupe, then skip it
                continue;
            } else if gp.contig == DUPE_NORMAL_LEVEL {
                let dupe_at_gp_pos = self.m_dupe_list.get(gp.position as usize).unwrap();
                for g in (0..dupe_at_gp_pos.len()) {
                    let gplong = gp_to_i64(&shift(dupe_at_gp_pos.get(g).unwrap(), i as i32));
                    // log::debug!("gplong={} else if", gplong);
                    kmer_stat.entry(gplong).and_modify(|s| *s += 1).or_insert(1);
                }
            } else {
                // not a dupe
                let gplong = gp_to_i64(&shift(gp, i as i32));
                // log::debug!("gplong={} else", gplong);
                kmer_stat.entry(gplong).and_modify(|s| *s += 1).or_insert(1);
            }
        }

        // get 1st and 2nd hit
        let mut gp1 = 0_i64;
        let mut count1 = 0;
        let mut gp2 = 0_i64;
        let mut count2 = 0;

        // log::debug!("kmer_stat={:#?}", kmer_stat);
        //TODO: handle small difference caused by INDEL
        // ocw.as_mut().and_then(|mut f| {
        //     writeln!(&mut f, "kmer_stat={:?}", kmer_stat).unwrap();
        //     Some(())
        // });

        for (k, v) in kmer_stat.iter() {
            if *k != 0 && *v > count1 {
                gp2 = gp1;
                count2 = count1;
                gp1 = *k;
                count1 = *v;
            } else if *k != 0 && *v > count2 {
                gp2 = *k;
                count2 = *v;
            }
        }

        // ocw.as_mut().and_then(|mut f| {
        //     writeln!(&mut f, "gp1={}, gp2={}, count1={}, count2={}", gp1, gp2, count1, count2).unwrap();
        //     Some(())
        // });

        if (count1 * step as i32) < global_settings().major_gene_key_requirement
            || (count2 * step as i32) < global_settings().minor_gene_key_requirement
        {
            // log::debug!("map_read, return empty vec. {}:{}", file!(), line!());
            // return an null list

            return Vec::new();
        }

        let mut mask: Vec<u8> = vec![MATCH_UNKNOWN; seqlen as usize];

        // second pass, make the mask
        kmer = -1;
        for i in (0..(seqlen as i32 - KMER + 1)) {
            // ocw.as_mut().and_then(|mut f| {
            //     write!(&mut f, "\"{}\",{},\"", "sp",i).unwrap();
            //     write!(&mut f, "\"\n", ).unwrap();
            //     // write_tsv_row!(&mut f, i, gp.contig, "2nd-pass");
            //     Some(())
            // });

            kmer = make_kmer_cv(&seq_cv, i, kmer, 1);
            // ocw.as_mut().and_then(|mut f| {
            //     write_tsv_row!(&mut f, i,kmer, "2nd-pass");
            //     Some(())
            // });
            if kmer < 0 {
                // ocw.as_mut().and_then(|mut f| {
                //     write!(&mut f, "\"k{}\",{},\"", kmer,i).unwrap();
                //     write!(&mut f, "\"\n", ).unwrap();
                //     // write_tsv_row!(&mut f, i, gp.contig, "2nd-pass");
                //     Some(())
                // });
                continue;
            }
            let pos = kmer.wrapping_shr(3);
            let bit = kmer & 0x07;

            // ocw.as_mut().and_then(|mut f| {
            //     log::debug!("i={}, bloom_filter_pos={}, bit={}", i, self.m_bloom_filter.get(pos as usize).unwrap(), 0x1_u8.wrapping_shl(bit as u32));
            //     Some(())
            // });

            if self.m_bloom_filter.get(pos as usize).unwrap()
                & (0x1_u8.wrapping_shl(bit as u32) as u8)
                == 0
            {
                // ocw.as_mut().and_then(|mut f| {
                //     let c = format!("{}, {}", self.m_bloom_filter.get(pos as usize).unwrap(), (0x1_u8.wrapping_shl(bit as u32) as u8));
                //     write!(&mut f, "\"b{}\",{},\"", c,i).unwrap();
                //     write!(&mut f, "\"\n", ).unwrap();
                //     // write_tsv_row!(&mut f, i, gp.contig, "2nd-pass");
                //     Some(())
                // });
                continue;
            }

            let gp = self.m_kmer_pos.get(&kmer).unwrap();
            // ocw.as_mut().and_then(|mut f| {
            //     write_tsv_row!(&mut f, i, gp.contig, "2nd-pass");
            //     Some(())
            // });
            // is a dupe
            if gp.contig == DUPE_HIGH_LEVEL {
                // ocw.as_mut().and_then(|mut f| {
                //     write!(&mut f, "\"c{}\",{},\"", gp.contig,i).unwrap();
                //     write!(&mut f, "\"\n", ).unwrap();
                //     // write_tsv_row!(&mut f, i, gp.contig, "2nd-pass");
                //     Some(())
                // });
                // too much keys in this dupe, then skip it
                continue;
            } else if gp.contig == DUPE_NORMAL_LEVEL {
                for g in (0..self.m_dupe_list.get(gp.position as usize).unwrap().len()) {
                    let gplong = gp_to_i64(&shift(
                        self.m_dupe_list
                            .get(gp.position as usize)
                            .unwrap()
                            .get(g)
                            .unwrap(),
                        i as i32,
                    ));

                    // ocw.as_mut().and_then(|mut f| {
                    //     write!(&mut f, "{},{},\"", g,i).unwrap();
                    //     write!(&mut f, "gplong={}\"\n", gplong).unwrap();
                    //     // write_tsv_row!(&mut f, i, gp.contig, "2nd-pass");
                    //     Some(())
                    // });

                    if (gplong - gp1).abs() <= 1 {
                        make_mask(mask.as_mut_slice(), MATCH_TOP, seqlen, i, KMER);
                        // ocw.as_mut().and_then(|mut f| {
                        //     write!(&mut f, "{},{},\"", g,i).unwrap();
                        //     for m in mask.iter() {
                        //         write!(&mut f, "{}", m).unwrap();
                        //     }
                        //     write!(&mut f, "\"\n", ).unwrap();
                        //     // write_tsv_row!(&mut f, i, gp.contig, "2nd-pass");
                        //     Some(())
                        // });
                    } else if (gplong - gp2).abs() <= 1 {
                        make_mask(mask.as_mut_slice(), MATCH_SECOND, seqlen, i, KMER);
                        // ocw.as_mut().and_then(|mut f| {
                        //     write!(&mut f, "{},{},\"", g,i).unwrap();
                        //     for m in mask.iter() {
                        //         write!(&mut f, "{}", m).unwrap();
                        //     }
                        //     write!(&mut f, "\"\n", ).unwrap();
                        //     // write_tsv_row!(&mut f, i, gp.contig, "2nd-pass");
                        //     Some(())
                        // });
                    } else if (gplong == 0) {
                        make_mask(mask.as_mut_slice(), MATCH_NONE, seqlen, i, KMER);
                        // ocw.as_mut().and_then(|mut f| {
                        //     write!(&mut f, "{},{},\"", g,i).unwrap();
                        //     for m in mask.iter() {
                        //         write!(&mut f, "{}", m).unwrap();
                        //     }
                        //     write!(&mut f, "\"\n", ).unwrap();
                        //     // write_tsv_row!(&mut f, i, gp.contig, "2nd-pass");
                        //     Some(())
                        // });
                    }
                }
            } else {
                let gplong = gp_to_i64(&shift(gp, i as i32));
                // ocw.as_mut().and_then(|mut f| {
                //     write!(&mut f, "{},{},\"", -1,i).unwrap();
                //     write!(&mut f, "gplong={}\"\n", gplong).unwrap();
                //     // write_tsv_row!(&mut f, i, gp.contig, "2nd-pass");
                //     Some(())
                // });
                if (gplong - gp1).abs() <= 1 {
                    make_mask(mask.as_mut_slice(), MATCH_TOP, seqlen, i, KMER);
                    // ocw.as_mut().and_then(|mut f| {
                    //     write!(&mut f, "-1,{},\"", i).unwrap();
                    //     for m in mask.iter() {
                    //         write!(&mut f, "{}", m).unwrap();
                    //     }
                    //     write!(&mut f, "\"\n", ).unwrap();
                    //     // write_tsv_row!(&mut f, i, gp.contig, "2nd-pass");
                    //     Some(())
                    // });
                } else if (gplong - gp2).abs() <= 1 {
                    make_mask(mask.as_mut_slice(), MATCH_SECOND, seqlen, i, KMER);
                    // ocw.as_mut().and_then(|mut f| {
                    //     write!(&mut f, "-1,{},\"", i).unwrap();
                    //     for m in mask.iter() {
                    //         write!(&mut f, "{}", m).unwrap();
                    //     }
                    //     write!(&mut f, "\"\n", ).unwrap();
                    //     // write_tsv_row!(&mut f, i, gp.contig, "2nd-pass");
                    //     Some(())
                    // });
                } else if gplong == 0 {
                    make_mask(mask.as_mut_slice(), MATCH_NONE, seqlen, i, KMER);
                    // ocw.as_mut().and_then(|mut f| {
                    //     write!(&mut f, "-1,{},\"", i).unwrap();
                    //     for m in mask.iter() {
                    //         write!(&mut f, "{}", m).unwrap();
                    //     }
                    //     write!(&mut f, "\"\n", ).unwrap();
                    //     // write_tsv_row!(&mut f, i, gp.contig, "2nd-pass");
                    //     Some(())
                    // });
                }
            }
        }

        let mut mismatches = 0;
        mask.iter_mut().take(seqlen as usize).for_each(|m| {
            if *m == MATCH_NONE || *m == MATCH_UNKNOWN {
                mismatches += 1;
            }
        });

        if mismatches > global_settings().mismatch_threshold {
            // too many mismatch indicates not a real fusion

            // log::debug!("map_read, return empty vec. {}:{}", file!(), line!());
            return Vec::new();
        }

        segment_mask(mask.as_slice(), seqlen, i64_to_gp(gp1), i64_to_gp(gp2))
    }

    /// this function is to gurantee that all the supporting reads will have same direction
    pub(crate) fn in_required_direction(&self, mapping: &[SeqMatch]) -> bool {
        if mapping.len() < 2 {
            return false;
        }

        let mut left = mapping.get(0).unwrap();
        let mut right = mapping.get(1).unwrap();

        if left.seq_start > right.seq_start {
            (left, right) = (right, left);
        }

        // both are positive, good to go
        if left.start_gp.position > 0 && right.start_gp.position > 0 {
            return true;
        }

        // if both are negative, we should use their reverse complement, which will be both positive
        if left.start_gp.position < 0 && right.start_gp.position < 0 {
            return false;
        }

        // if one is positive, the other is negative, their reverse complement will be the same
        if self
            .m_fusions
            .get(left.start_gp.contig as usize)
            .unwrap()
            .is_reversed()
            && !self
                .m_fusions
                .get(right.start_gp.contig as usize)
                .unwrap()
                .is_reversed()
        {
            // if left is reversed gene and right is forward gene
            // we should use their reverse complement, which keep the left forward gene
            return false;
        } else if !self
            .m_fusions
            .get(left.start_gp.contig as usize)
            .unwrap()
            .is_reversed()
            && self
                .m_fusions
                .get(right.start_gp.contig as usize)
                .unwrap()
                .is_reversed()
        {
            // if left is forward gene and right is reversed gene, good to go
            return true;
        } else {
            // otherwise, we should keep the left has smaller contig
            if left.start_gp.contig < right.start_gp.contig {
                return true;
            }
            // or smaller positive if contig is the same
            if left.start_gp.contig == right.start_gp.contig
                && left.start_gp.position.abs() < left.start_gp.position.abs()
            //TODO: left  < left ?
            {
                return true;
            } else {
                return false;
            }
        }

        false
    }

    fn print_stat(&self) {
        println!("m_unique_pos:{}", self.m_unique_pos);
        println!("m_dupe_pos:{}", self.m_dupe_pos);
    }
}

fn segment_mask(mask: &[u8], seqlen: i32, gp1: GenePos, gp2: GenePos) -> Vec<SeqMatch> {
    let mut result: Vec<SeqMatch> = Vec::new();

    const ALLOWED_GAP: i32 = 10;
    const THRESHOLD_LEN: i32 = 20;

    let targets = [MATCH_TOP as i32, MATCH_SECOND as i32];
    let gps = [gp1, gp2];

    for (target, gp) in targets.into_iter().zip(gps.into_iter()) {
        let mut max_start = -1_i32;
        let mut max_end = -1_i32;

        // get gp1
        let mut start = 0_i32;
        let mut end = 0_i32;

        loop {
            // get next start
            while *mask.get(start as usize).unwrap() as i32 != target && start != seqlen - 1 {
                start += 1;
            }
            // reach the tail
            if start >= seqlen - 1 {
                break;
            }

            if *mask.get(start as usize).unwrap() as i32 == target {
                end = start + 1;
                // get the end
                let mut g = 0_i32;
                while (g < ALLOWED_GAP) && (end + g) < seqlen {
                    if *mask.get((end + g) as usize).unwrap() as i32 > target {
                        break;
                    }

                    if end + g < seqlen && *mask.get((end + g) as usize).unwrap() as i32 == target {
                        end += g + 1;
                        g = 0;
                        continue;
                    }

                    g += 1;
                }
                // left shift to remove the mismatched end
                end -= 1;
                if end - start > (max_end - max_start) {
                    max_end = end;
                    max_start = start;
                }
                start += 1;
            } else {
                // not found
                break;
            }
        }
        if max_end - max_start > THRESHOLD_LEN {
            let seq_match = SeqMatch::new(max_start, max_end, gp);
            result.push(seq_match);
        }
    }

    result
}

#[inline]
fn concat_i32_bits_into_i64(two_i32_array: [i32; 2]) -> i64 {
    let mut arr_iter = two_i32_array.into_iter().flat_map(|n| n.to_ne_bytes());
    let concated_ne_bytes: [u8; 8] = array::from_fn(|_| arr_iter.next().unwrap());

    i64::from_ne_bytes(concated_ne_bytes)
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

    let concated_bits = concat_i32_bits_into_i64([gp.position, 0]);

    (ret.wrapping_shl(32)) | concated_bits

    // gp.position
}

#[inline]
pub(crate) fn i64_to_gp(val: i64) -> GenePos {
    GenePos {
        contig: (val.wrapping_shr(32)) as i16,
        position: (val & 0x00000000FFFFFFFF) as i32,
    }
}

fn make_mask(mask: &mut [u8], flag: u8, seqlen: i32, start: i32, kmer_size: i32) {
    let end_point = seqlen.min((start + kmer_size));

    // mask.iter_mut()
    //     .skip(start as usize)
    //     .take((end_point - start) as usize)
    //     .for_each(|m| {
    //         *m = (*m).max(flag);
    //     })

    mask.get_mut((start as usize)..(end_point as usize))
        .unwrap()
        .iter_mut()
        .for_each(|m| {
            *m = (*m).max(flag);
        })
}

fn make_kmer(seq: &str, pos: i32, last_kmer: i32, step: i32) -> i32 {
    // else calculate it completely
    let mut kmer = 0;

    let mut start = 0;

    // re-use several bits
    if last_kmer >= 0 {
        kmer = last_kmer;
        start = KMER - step;

        if step == 1 {
            kmer = (kmer & 0x3FFFFFFF).wrapping_shl(2);
        } else if step == 2 {
            kmer = (kmer & 0x0FFFFFFF).wrapping_shl(2);
        } else if step == 3 {
            kmer = (kmer & 0x03FFFFFF).wrapping_shl(2);
        } else if step == 4 {
            kmer = (kmer & 0x00FFFFFF).wrapping_shl(2);
        }
    }

    for (base, i) in seq
        .chars()
        .skip((pos + start) as usize)
        .take((KMER - start) as usize)
        .zip(start..KMER)
    {
        match base {
            'A' => {
                kmer += 0;
            }
            'T' => {
                kmer += 1;
            }
            'C' => {
                kmer += 2;
            }
            'G' => {
                kmer += 3;
            }
            _ => {
                return -1;
            }
        }

        // not the tail
        if (i < KMER - 1) {
            kmer = kmer << 2;
        }
    }

    kmer
}

fn make_kmer_cv(seq: &[char], pos: i32, last_kmer: i64, step: i32) -> i64 {
    // else calculate it completely
    let mut kmer = 0;

    let mut start = 0;

    // re-use several bits
    if last_kmer >= 0 {
        kmer = last_kmer;
        start = KMER - step;

        if step == 1 {
            kmer = (kmer & 0x3FFFFFFF).wrapping_shl(2);
        } else if step == 2 {
            kmer = (kmer & 0x0FFFFFFF).wrapping_shl(2);
        } else if step == 3 {
            kmer = (kmer & 0x03FFFFFF).wrapping_shl(2);
        } else if step == 4 {
            kmer = (kmer & 0x00FFFFFF).wrapping_shl(2);
        }
    }
    let subseq = seq
        .get(((pos + start) as usize)..((pos + KMER) as usize))
        .unwrap_or_else(|| {
            panic!(
                "seq_cv_len={}, pos={}, s={}, e={}",
                seq.len(),
                pos,
                pos + start,
                pos + KMER
            )
        });

    // log::debug!("subseq.len() = {}", subseq.len());
    // log::debug!("pos={}, start={}, KMER = {}", pos, start, KMER);
    for (base, i) in subseq.iter().zip(start..KMER) {
        match base {
            'A' => {
                kmer += 0;
            }
            'T' => {
                kmer += 1;
            }
            'C' => {
                kmer += 2;
            }
            'G' => {
                kmer += 3;
            }
            _ => {
                return -1;
            }
        }

        // not the tail
        if (i < KMER - 1) {
            kmer = kmer << 2;
        }
    }

    kmer
}

fn make_kmer_bytes(seq: &str, pos: i32, last_kmer: i64, step: i32) -> i64 {
    // else calculate it completely
    let mut kmer = 0;

    let mut start = 0;

    // re-use several bits
    if last_kmer >= 0 {
        kmer = last_kmer;
        start = KMER - step;

        if step == 1 {
            kmer = (kmer & 0x3FFFFFFF).wrapping_shl(2);
        } else if step == 2 {
            kmer = (kmer & 0x0FFFFFFF).wrapping_shl(2);
        } else if step == 3 {
            kmer = (kmer & 0x03FFFFFF).wrapping_shl(2);
        } else if step == 4 {
            kmer = (kmer & 0x00FFFFFF).wrapping_shl(2);
        }
    }
    let subseq = seq
        .get(((pos + start) as usize)..((pos + KMER) as usize))
        .unwrap_or_else(|| {
            panic!(
                "seq_cv_len={}, pos={}, s={}, e={}",
                seq.len(),
                pos,
                pos + start,
                pos + KMER
            )
        });

    // log::debug!("subseq.len() = {}", subseq.len());
    // log::debug!("pos={}, start={}, KMER = {}", pos, start, KMER);
    for (base, i) in subseq.as_bytes().iter().zip(start..KMER) {
        match base {
            b'A' => {
                kmer += 0;
            }
            b'T' => {
                kmer += 1;
            }
            b'C' => {
                kmer += 2;
            }
            b'G' => {
                kmer += 3;
            }
            _ => {
                return -1;
            }
        }

        // not the tail
        if (i < KMER - 1) {
            kmer = kmer << 2;
        }
    }

    kmer
}

#[inline(always)]
fn check_map_has_key_then_get_back<K, V>(map: &BTreeMap<K, V>, key: K) -> Option<K>
where
    K: Hash + std::cmp::Eq + Ord,
{
    match map.contains_key(&key) {
        true => Some(key),
        false => None,
    }
}

#[cfg(test)]
mod test {
    use std::{array, ops::Shl};

    use crate::{
        core::{
            common::GenePos, fusion::Fusion, indexer::concat_i32_bits_into_i64, read::SequenceRead,
            sequence::Sequence,
        },
        utils::logging::init_logger,
    };

    use super::{gp_to_i64, i64_to_gp, make_kmer_cv, Indexer};

    #[test]
    fn bit1() {
        println!("{:0>32b}", 10_i32,);
        println!("{:b}", -10_i32,);

        // let c = -1_i32.shl(32_i32);
        // println!("{:0>32b} {}", c, c);

        let c = 20_i32.wrapping_shl(32);
        println!("{:0>32b} {}", c, c);

        // let c = [1231231, 0] as i64;
        // println!("{:0>64b} {}", c, c);
    }

    #[test]
    fn concat_bits() {
        let a = [1234351_i32, 0];

        for e in a.iter() {
            println!("{:0>32b}", e);
        }

        // let mut iter1 = a.iter().flat_map(|e| e.to_ne_bytes());

        // let arr1: [u8; 8] = array::from_fn(|_| iter1.next().unwrap());

        // let b = i64::from_ne_bytes(arr1);

        let b = concat_i32_bits_into_i64(a);

        println!("{:0>64b}", b);

        println!("{:0>64b}", b & 0x00000000FFFFFFFF);
    }

    #[test]
    fn convert() {
        println!("{}", _convert());
    }

    fn _convert() -> bool {
        let contigs = [0, 1, 3, 220, -1, 0, 23, 4440, 110, 10];
        let positions = [0, 111, 222, -333, 444, 555555, 6, -7777777, 8888, -9999];

        for i in (0..10) {
            let gp = GenePos {
                contig: contigs[i],
                position: positions[i],
            };
            let val = gp_to_i64(&gp);
            println!("val={:0>64b}", val);
            let gp2 = i64_to_gp(val);

            let val2 = gp_to_i64(&gp2);
            println!("val2={:0>64b}", val2);

            if (gp.contig == gp2.contig && gp.position == gp2.position && val == val2) == false {
                println!(
                    "{:#?}",
                    (
                        i,
                        gp.contig,
                        gp2.contig,
                        gp.position,
                        gp2.position,
                        val,
                        val2
                    )
                );

                return false;
            }
        }

        true
    }

    #[test]
    fn mut_ref() {
        let mut a = 1;
        let ar = &mut a;
        let b = (*ar).max(3);

        *ar = b;

        println!("{}", b);
    }

    #[test]
    fn test_indexing() {
        init_logger();
        let fusion_list = Fusion::parse_csv(
            "/home/eck/workspace/240115_genefuse_rust/genefuse_rust/genes/druggable.hg38.csv",
        )
        .unwrap();
        let mut m_indexer = Indexer::new(
            "/home/eck/workspace/240115_genefuse_rust/genefuse_rust/hg38.fa.gz",
            fusion_list.clone(),
        )
        .unwrap();
        m_indexer.make_index();
    }

    #[test]
    fn indexer_map_read() {
        init_logger();
        let fusion_list = Fusion::parse_csv(
            "/home/eck/workspace/240115_genefuse_rust/genefuse_rust/genes/druggable.hg38.csv",
        )
        .unwrap();
        let mut m_indexer = Indexer::new(
            "/home/eck/workspace/240115_genefuse_rust/genefuse_rust/hg38.fa.gz",
            fusion_list.clone(),
        )
        .unwrap();

        m_indexer.make_index();

        let s = SequenceRead { m_name: "@NB551106:23:HVMTYBGX2:2:12302:19642:13894 1:N:0:GATCAG merged_diff_0".into(), m_seq: Sequence { m_str: "CATCACACACCTTGACTGGTCCCCAGACAACAAGTATATAATGTCTAACTCGGGAGACTATGAAATATTGTACTGTAAGTATGAATGATTTTATATATATATATATATGCTATGATTATATTTATATATATAATAATTATTTTCCATATATCTGATTTTTAGCTTTGCATTTACTTTA".into() }, m_strand: "+".into(), m_quality: "A/AA/EEEEEAEEAEEEEEEAEEEAEAZZSZZZDZZZZZZZZZZZZZZZZZZSZZZZZZZZZZZZSZSSOZZZZZSSZZZZZZZZZZZVSSSZZZZZSZZZZZZZZZZSZZZZZZSZSZZZ=SZZZZZZZ=ZSZZZZSZZ=SZZZZJZZZZEEEA66AEEEEE666AEE6EEEA6A/A".into(), m_has_quality: true };

        println!("{:?}", m_indexer.map_read(&s));
    }

    #[test]
    fn use_parenthesis() {
        println!("{}", 254_i32 - 234_u32 as i32 == 20_i32);
        println!("{}", (254_i32 - 234_u32 as i32) == 20_i32);
    }

    #[test]
    fn make_kmer_test() {
        init_logger();
        let seq = "CATCACACACCTTGACTGGTCCCCAGACAACAAGTATATAAT"
            .chars()
            .collect::<Vec<_>>();

        println!("{}", seq.get(0..16).unwrap().len());
        println!("{}", make_kmer_cv(&seq, 0, -1, 2));
    }

    #[test]
    fn bool_cmp() {
        println!("{}", 96 & 64 as u8);
        println!("{}", 96 & 64);
        println!("{}", 96 & 64 as u8 == 0);
    }
}
