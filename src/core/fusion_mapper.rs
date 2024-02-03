use std::{
    borrow::Cow,
    cmp::{Ordering, Reverse},
    error::Error,
    sync::Mutex,
};

use crate::{
    aux::{
        global_settings::global_settings,
        pbar::{prepare_pbar, PBSummary},
    },
    core::{fusion, pescanner::DBT, sequence::Sequence},
    utils::{dis_connected_count, StringCPP},
};

use super::{
    edit_distance::edit_distance,
    fasta_reader::FastaReader,
    fusion::Fusion,
    fusion_result::FusionResult,
    indexer::{Indexer, SeqMatch},
    matcher::Matcher,
    read::SequenceRead,
    read_match::ReadMatch,
    sequence::reverse_complement,
};

pub(crate) struct FusionMapper {
    pub(crate) m_ref_file: String,
    pub(crate) m_fusion_match_size: i32,
    pub(crate) m_indexer: Indexer,
    pub(crate) fusion_list: Vec<Fusion>,
    pub(crate) fusion_matches: Mutex<Vec<Vec<ReadMatch>>>,
    pub(crate) m_fusion_results: Vec<FusionResult>,
}

impl FusionMapper {
    pub(crate) fn from_ref_and_fusion_files(
        ref_file: &str,
        fusion_file: &str,
    ) -> Result<Self, Box<dyn Error>> {
        let fusion_list = Fusion::parse_csv(fusion_file)?;
        log::debug!("Parsed csv.");

        // log::debug!("fusion_list={:#?}", fusion_list);

        let mut m_indexer = Indexer::new(ref_file, fusion_list.clone())?;
        m_indexer.make_index();
        log::debug!("Made index.");

        // init()
        let m_fusion_match_size = fusion_list.len().pow(2);

        let fusion_matches = Mutex::new(vec![vec![]; m_fusion_match_size]);

        Ok(Self {
            m_ref_file: ref_file.to_string(),
            m_fusion_match_size: m_fusion_match_size as i32,
            m_indexer,
            fusion_list,
            fusion_matches,
            m_fusion_results: Vec::new(),
        })
    }

    pub(crate) fn from_fasta_reader_and_fusion_files(
        fasta_reader:FastaReader,
        fusion_file: &str,
    ) -> Result<Self, Box<dyn Error>> {
        let fusion_list = Fusion::parse_csv(fusion_file)?;
        log::debug!("Parsed csv.");

        // log::debug!("fusion_list={:#?}", fusion_list);

        let mut m_indexer = Indexer::with_loaded_ref(fasta_reader, fusion_list.clone());
        m_indexer.make_index();
        log::debug!("Made index.");

        // init()
        let m_fusion_match_size = fusion_list.len().pow(2);

        let fusion_matches = Mutex::new(vec![vec![]; m_fusion_match_size]);

        Ok(Self {
            m_ref_file: m_indexer.m_reference.as_ref().unwrap().m_fasta_file.clone(),
            m_fusion_match_size: m_fusion_match_size as i32,
            m_indexer,
            fusion_list,
            fusion_matches,
            m_fusion_results: Vec::new(),
        })
    }

    /// ## Defaults:
    /// - distance_req = 2
    ///
    /// - qual_req = 20
    pub(crate) fn map_read(
        &self,
        r: &SequenceRead,
        mapable: &mut bool,
        distance_req: i32,
        qual_req: i32,
    ) -> Result<Option<ReadMatch>, Box<dyn Error>> {
        let mut mapping = self.m_indexer.map_read(r);

        // if r.m_name.contains(DBT) {
        //     log::debug!("mapping={:#?}", mapping);
        // }

        //we only focus on the reads that can be mapped to two genome positions
        if mapping.len() < 2 {
            // if r.m_name.contains(DBT) {
            //     log::debug!("mapping.len()={}", mapping.len());
            // }
            *mapable = false;
            return Ok(None);
        }

        *mapable = true;

        //if the left part of mapping result is reverse, use its reverse complement alternative and skip this one
        if !self.m_indexer.in_required_direction(&mapping) {
            // if r.m_name.contains(DBT) {
            //     log::debug!("in_required_direction = false");
            // }
            return Ok(None);
        }

        // TODO: set int readBreak, int leftContig, int leftPos, int rightContig, int rightPos
        let m = self.make_match(r, &mut mapping);
        // if r.m_name.contains(DBT) {
        //     log::debug!("make_match res = {:#?}", m);
        // }

        Ok(m)
    }

    fn get_ref(&self) -> Option<&FastaReader> {
        // indexer can be NULL in cpp code.
        // if self.m_indexer.is_none() {
        //     None
        // } else {
        //     self.m_indexer.get_ref()
        // }
        self.m_indexer.get_ref()
    }

    fn get_ref_mut(&mut self) -> Option<&mut FastaReader> {
        // indexer can be NULL in cpp code.
        // if self.m_indexer.is_none() {
        //     None
        // } else {
        //     self.m_indexer.get_ref()
        // }
        self.m_indexer.get_ref_mut()
    }

    fn make_match(&self, r: &SequenceRead, mapping: &mut [SeqMatch]) -> Option<ReadMatch> {
        if mapping.len() != 2 {
            return None;
        }
        let mapping_split = mapping.split_first_mut().unwrap();

        let mut left = mapping_split.0;
        let mut right = mapping_split.1.first_mut().unwrap();

        if left.seq_start > right.seq_start {
            (left, right) = (right, left);
        }

        // log::debug!(
        //     "left.seq_start={}, right.seq_start={}",
        //     left.seq_start,
        //     right.seq_start
        // );

        let read_break = (left.seq_end + right.seq_start) / 2;
        let left_gp = &mut left.start_gp;
        let right_gp = &mut right.start_gp;

        left_gp.position += read_break;
        right_gp.position += read_break + 1;

        let gap = right.seq_start - left.seq_end - 1;

        let mut read_match = ReadMatch::new(
            r.clone(),
            read_break,
            left_gp.clone(),
            right_gp.clone(),
            gap,
            false,
        );

        self.calc_distance(&mut read_match);

        Some(read_match)
    }

    fn calc_distance(&self, m: &mut ReadMatch) -> () {
        let seq = m.m_read.m_seq.m_str.as_str();

        let read_break = m.m_read_break;
        let left_len = read_break + 1;
        let right_len = seq.len() as i32 - (read_break + 1);

        let left_seq = seq.subchars(0, left_len as usize);
        let right_seq = seq.subchars((read_break + 1) as usize, right_len as usize);

        // let left_gene = &self.fusion_list.get(m.m_left_gp.contig as usize).unwrap().m_gene;
        // let right_gene = &self.fusion_list.get(m.m_right_gp.contig as usize).unwrap().m_gene;

        m.m_left_distance = self.calc_ed(
            &left_seq,
            m.m_left_gp.contig as i32,
            m.m_left_gp.position - left_len + 1,
            m.m_left_gp.position,
        );

        m.m_right_distance = self.calc_ed(
            &right_seq,
            m.m_right_gp.contig as i32,
            m.m_right_gp.position,
            m.m_right_gp.position + right_len - 1,
        );
    }

    fn calc_ed(&self, seq: &str, contig: i32, mut start: i32, mut end: i32) -> i32 {
        // check start and end are in same strand
        if (start >= 0 && end <= 0) || (start <= 0 && end >= 0) {
            return -1;
        }

        let fusion_seq = self.m_indexer.m_fusion_seq.get(contig as usize).unwrap();

        // check the overflow
        if start.abs() >= (fusion_seq.len() as i32) || end.abs() >= (fusion_seq.len() as i32) {
            return -2;
        }

        let mut ss = Cow::from(seq);
        if start < 0 {
            let s = Sequence::new(seq.to_owned());
            let rc = s.reverse_complement();
            ss = Cow::from(reverse_complement(seq));

            let tmp = start;
            start = -end;
            end = -tmp;
        }

        let ref_str = fusion_seq.subchars(start as usize, (end - start + 1) as usize);

        edit_distance(&ss, ss.len(), &ref_str, ref_str.len()) as i32
    }

    pub(crate) fn add_match(&self, m: ReadMatch) -> () {
        let left_contig = m.m_left_gp.contig;
        let right_contig = m.m_right_gp.contig;

        log::debug!(
            "add_match(), left_contig={}, right_contig={}",
            left_contig,
            right_contig
        );

        let index = self.fusion_list.len() as i32 * right_contig as i32 + left_contig as i32;

        let mut fusion_matches = self.fusion_matches.lock().unwrap();

        match fusion_matches.get_mut(index as usize) {
            Some(v) => v.push(m),
            None => panic!(
                "index={}, fusion_mathces_len={}",
                index,
                fusion_matches.len()
            ),
        }
    }
    pub(crate) fn filter_matches(&mut self) -> () {
        // calc the sequence number before any filtering
        // let mut total = 0;
        // for fm in self.fusion_matches.lock().unwrap().iter() {
        //     total += fm.len();
        // }

        let total = self
            .fusion_matches
            .lock()
            .unwrap()
            .iter()
            .fold(0, |a, b| a + b.len());

        log::info!("sequence number before filtering: {}", total);

        self.remove_by_complexity();
        self.remove_by_distance();
        self.remove_indels();
        self.remove_alignables();
    }

    fn remove_by_complexity(&mut self) -> () {
        let mut removed = 0;
        for fm in self.fusion_matches.lock().unwrap().iter_mut() {
            fm.retain(|rm| {
                let seq = rm.m_read.m_seq.m_str.as_str();
                let read_break = rm.m_read_break;

                let dec = is_low_complexity(&seq.subchars(0, (read_break + 1) as usize))
                    || is_low_complexity(&seq.subchars(
                        (read_break + 1) as usize,
                        seq.len() - (read_break + 1) as usize,
                    ));

                if dec {
                    removed += 1;
                    false
                } else {
                    true
                }
            })
        }

        log::info!("remove_by_complexity: {}", removed);
    }

    fn remove_by_distance(&mut self) -> () {
        // diff should be less than DIFF_THRESHOLD
        const DIFF_THRESHOLD: i32 = 5;

        let mut removed = 0;
        self.fusion_matches
            .lock()
            .unwrap()
            .iter_mut()
            .for_each(|fm| {
                {
                    fm.retain(|rm| {
                        let dec = rm.m_left_distance + rm.m_right_distance >= DIFF_THRESHOLD;

                        if dec {
                            removed += 1;
                            false
                        } else {
                            true
                        }
                    })
                }
            });

        log::info!("removeByDistance: {}", removed);
    }

    fn remove_indels(&mut self) -> () {
        // diff should be greather than INDEL_THRESHOLD
        let INDEL_THRESHOLD = global_settings().deletion_threshold;
        let mut removed = 0;

        self.fusion_matches
            .lock()
            .unwrap()
            .iter_mut()
            .for_each(|fm| {
                {
                    fm.retain(|rm| {
                        let dec = rm.m_left_gp.contig == rm.m_right_gp.contig
                            && (rm.m_left_gp.position - rm.m_right_gp.position).abs()
                                < INDEL_THRESHOLD as i32;

                        if dec {
                            removed += 1;
                            false
                        } else {
                            true
                        }
                    })
                }
            });

        log::info!("removeIndels: {}", removed);
    }

    pub(crate) fn sort_matches(&mut self) {
        // sort the matches to make the pileup more clear
        self.fusion_matches
            .lock()
            .unwrap().iter_mut().for_each(|rmv|{
                rmv.sort_by(|a, b| b.partial_cmp(a).unwrap());
            });
            




        
    }

    pub(crate) fn free_matches(&mut self) {
        // free it
        self.fusion_matches.lock().unwrap().clear();
    }

    pub(crate) fn cluster_matches(&mut self) {
        log::debug!("self.m_fusion_match_size={}", self.m_fusion_match_size);
        log::debug!(
            "fusion_matches_len={}",
            self.fusion_matches.lock().unwrap().len()
        );
        for (i, fm) in (0..(self.m_fusion_match_size)).zip(
            self.fusion_matches
                .lock()
                .unwrap()
                .iter()
                .map(|e| e.as_slice()),
        ) {
            let mut frs = Vec::<FusionResult>::new();
            for (m, rm) in (0..fm.len()).zip(fm.iter()) {
                let mut found = false;
                for (f, fr) in (0..frs.len()).zip(frs.iter_mut()) {
                    if fr.support(rm) {
                        fr.add_match(rm.clone());
                        found = true;
                        break;
                    }
                }
                if !found {
                    let mut fr = FusionResult::with_minimum();
                    fr.add_match(rm.clone());
                    frs.push(fr);
                }
            }
            for (f, mut fr) in (0..frs.len()).zip(frs.into_iter()) {
                // log::debug!("init -> fr.m_left_ref_ext={} fr.m_right_ref={}", fr.m_left_ref_ext, fr.m_right_ref);
                log::debug!("init -> fusion_list={:#?}", self.fusion_list);
                fr.calc_fusion_point();
                log::debug!("calc_fusion_point -> fusion_list={:#?}", self.fusion_list);

                // log::debug!("calc_fusion_point -> fr.m_left_ref_ext={} fr.m_right_ref={}", fr.m_left_ref_ext, fr.m_right_ref);
                // log::debug!("fr.m_left_gp.contig={}, fr.m_right_gp.contig={}", fr.m_left_gp.contig, fr.m_right_gp.contig);
                fr.make_reference(
                    self.m_indexer
                        .m_fusion_seq
                        .get(fr.m_left_gp.contig as usize)
                        .unwrap_or_else(|| {
                            panic!(
                                "self.m_indexer.m_fusion_seq.len()={}, index={}",
                                self.m_indexer.m_fusion_seq.len(),
                                fr.m_left_gp.contig
                            )
                        }),
                    self.m_indexer
                        .m_fusion_seq
                        .get(fr.m_right_gp.contig as usize)
                        .unwrap(),
                );
                log::debug!("make_reference -> fusion_list={:#?}", self.fusion_list);
                // log::debug!("make_reference -> fr.m_left_ref_ext={} fr.m_right_ref={}", fr.m_left_ref_ext, fr.m_right_ref);
                fr.adjust_fusion_break();
                log::debug!("adjust_fusion_break -> fusion_list={:#?}", self.fusion_list);
                // log::debug!("adjust_fusion_break -> fr.m_left_ref_ext={} fr.m_right_ref={}", fr.m_left_ref_ext, fr.m_right_ref);
                fr.calc_unique();
                log::debug!("calc_unique -> fusion_list={:#?}", self.fusion_list);
                // log::debug!("calc_unique -> fr.m_left_ref_ext={} fr.m_right_ref={}", fr.m_left_ref_ext, fr.m_right_ref);
                fr.update_info(&self.fusion_list);
                log::debug!("update_info -> fusion_list={:#?}", self.fusion_list);
                // log::debug!("update_info -> fr.m_left_ref_ext={} fr.m_right_ref={}", fr.m_left_ref_ext, fr.m_right_ref);
                if fr.is_qualified() {
                    if !global_settings().output_deletions && fr.is_deletion() {
                        continue;
                    }
                    if fr.is_left_protein_forward() != fr.is_right_protein_forward() {
                        if !global_settings().output_untranslated {
                            continue;
                        }
                    }
                    fr.print(&self.fusion_list);
                    self.m_fusion_results.push(fr);
                }
            }
        }

        self.sort_fusion_results();
        log::info!("found {} fusions", self.m_fusion_results.len(),);
    }

    fn remove_alignables(&mut self) -> () {
        if self.get_ref().is_none() {
            return;
        }

        // let fasta_ref = fasta_refo.unwrap();

        let mut seqs: Vec<Sequence> = Vec::new();

        // first pass to gather all sequences
        for fm in self.fusion_matches.lock().unwrap().iter_mut() {
            for rm in fm.iter() {
                seqs.push(rm.get_read().m_seq.clone());
            }
        }

        log::info!("making matcher...");
        let mut matcher = Matcher::from_ref_and_seqs(self.m_indexer.get_ref_mut(), &seqs);

        let mut removed = 0;

        log::info!("removing alignable sequences...");
        // second pass to remove alignable sequences
        {
            let mut fusion_matches = self.fusion_matches.lock().unwrap();
            // let pb = prepare_pbar(fusion_matches.len() as u64);
            // pb.set_message("do_matching...");

            fusion_matches.iter_mut().for_each(|fm| {
                {
                    // let pb = prepare_pbar(fm.len() as u64);
                    // pb.set_message("do_matching for a fusion match ...");
                    
                    fm.retain(|rm| {
                        let mr = matcher.do_match(&rm.get_read().m_seq);
                        let dec = mr.is_some();

                        // pb.inc(1);
                        if dec {
                            removed += 1;
                            false
                        } else {
                            true
                        }
                    });

                    // pb.inc(1);
                    // pb.finish_and_clear();
                }
            });

            // pb.finish();
        }
        log::info!("removeAlignables: {}", removed);
    }

    fn sort_fusion_results(&mut self) {
        self.m_fusion_results
            .sort_by(|a, b| Self::more_reads(b, a).unwrap()) // b,a instead of a,b because we want descending order.
    }

    fn more_reads(r1: &FusionResult, r2: &FusionResult) -> Option<Ordering> {
        match r1.m_unique.partial_cmp(&r2.m_unique) {
            Some(Ordering::Equal) => {}
            ord => return ord,
        }

        r1.m_matches.len().partial_cmp(&r2.m_matches.len())
    }
}

fn is_low_complexity(s: &str) -> bool {
    if s.len() < 20 {
        return true;
    }

    if dis_connected_count(s) < 7 {
        return true;
    }

    false
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn sof() {
        // FusionMapper::from_ref_and_fusion_files(
        //     "/home/eck/workspace/240115_genefuse_rust/genefuse_rust/hg38.fa",
        //     "/home/eck/workspace/240115_genefuse_rust/genefuse_rust/genes/druggable.hg38.csv",
        // ).unwrap();
        // let a:Vec<Vec<ReadMatch>> = Vec::new();
        // let b = Indexer::new(
        //     "/home/eck/workspace/240115_genefuse_rust/genefuse_rust/hg38.fa",
        //     Vec::new(),
        // );
        const BLOOM_FILTER_SIZE: usize = (1_i32.wrapping_shl(20)) as usize;

        let c = [0_u8; BLOOM_FILTER_SIZE];
    }
    #[test]
    fn pcmp() {
        println!("{:?}", 1_i32.partial_cmp(&2).unwrap());

        let mut a = [5, 4,1,2,3];
        a.sort_by(|a,b| {
            let r = a.partial_cmp(b).unwrap();
            println!("{:?}, a={}, b={}", r, a,b);
            r
        });
        println!("a={:?}", a);
    }
}
