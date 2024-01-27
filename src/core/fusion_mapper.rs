use std::{
    borrow::Cow,
    cmp::{Ordering, Reverse},
    error::Error,
};

use crate::{
    aux::global_settings::global_settings,
    core::{fusion, sequence::Sequence},
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
    pub(crate) fusion_matches: Vec<Vec<ReadMatch>>,
    pub(crate) m_fusion_results: Vec<FusionResult>,
}

impl FusionMapper {
    pub(crate) fn from_ref_and_fusion_files(
        ref_file: &str,
        fusion_file: &str,
    ) -> Result<Self, Box<dyn Error>> {
        let fusion_list = Fusion::parse_csv(fusion_file)?;
        log::debug!("Parsed csv.");

        
        let mut m_indexer = Indexer::new(ref_file, fusion_list.clone())?;
        m_indexer.make_index();
        log::debug!("Made index.");

        let m_fusion_match_size = fusion_list.len().pow(2);

        let fusion_matches = Vec::with_capacity(m_fusion_match_size);

        Ok(Self {
            m_ref_file: ref_file.to_string(),
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
        &mut self,
        r: &SequenceRead,
        mapable: &mut bool,
        distance_req: i32,
        qual_req: i32,
    ) -> Result<Option<ReadMatch>, Box<dyn Error>> {
        let mut mapping = self.m_indexer.map_read(r);

        //we only focus on the reads that can be mapped to two genome positions
        if mapping.len() < 2 {
            *mapable = false;
            return Ok(None);
        }

        *mapable = true;

        //if the left part of mapping result is reverse, use its reverse complement alternative and skip this one
        if self.m_indexer.in_required_direction(&mapping) {
            return Ok(None);
        }

        // TODO: set int readBreak, int leftContig, int leftPos, int rightContig, int rightPos
        let m = self.make_match(r, &mut mapping);

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

    pub(crate) fn add_match(&mut self, m: ReadMatch) -> () {
        let left_config = m.m_left_gp.contig;
        let right_config = m.m_right_gp.contig;

        let index = self.fusion_list.len() as i32 * right_config as i32 + left_config as i32;

        self.fusion_matches.get_mut(index as usize).unwrap().push(m);
    }

    pub(crate) fn filter_matches(&mut self) -> () {
        // calc the sequence number before any filtering
        let mut total = 0;
        for fm in self.fusion_matches.iter() {
            total += fm.len();
        }

        log::info!("sequence number before filtering: {}", total);

        self.remove_by_complexity();
        self.remove_by_distance();
        self.remove_indels();
        self.remove_alignables();
    }

    fn remove_by_complexity(&mut self) -> () {
        let mut removed = 0;
        for fm in self.fusion_matches.iter_mut() {
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
        self.fusion_matches.iter_mut().for_each(|fm| {
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

        self.fusion_matches.iter_mut().for_each(|fm| {
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
            .sort_by(|a, b| b.partial_cmp(a).unwrap());
    }

    pub(crate) fn free_matches(&mut self) {
        // free it
        self.fusion_matches.clear();
    }

    pub(crate) fn cluster_matches(&mut self) {
        for (i, fm) in
            (0..(self.m_fusion_match_size)).zip(self.fusion_matches.iter().map(|e| e.as_slice()))
        {
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
                    let mut fr = FusionResult::default();
                    fr.add_match(rm.clone());
                    frs.push(fr);
                }
            }
            for (f, mut fr) in (0..frs.len()).zip(frs.into_iter()) {
                fr.calc_fusion_point();
                fr.make_reference(
                    self.m_indexer
                        .m_fusion_seq
                        .get(fr.m_left_gp.contig as usize)
                        .unwrap(),
                    self.m_indexer
                        .m_fusion_seq
                        .get(fr.m_right_gp.contig as usize)
                        .unwrap(),
                );
                fr.adjust_fusion_break();
                fr.calc_unique();
                fr.update_info(&self.fusion_list);
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
        for fm in self.fusion_matches.iter_mut() {
            for rm in fm.iter() {
                seqs.push(rm.get_read().m_seq.clone());
            }
        }

        let mut matcher = Matcher::from_ref_and_seqs(self.m_indexer.get_ref_mut(), &seqs);

        let mut removed = 0;
        // second pass to remove alignable sequences
        self.fusion_matches.iter_mut().for_each(|fm| {
            {
                fm.retain(|rm| {
                    let mr = matcher.do_match(&rm.get_read().m_seq);
                    let dec = mr.is_some();

                    if dec {
                        removed += 1;
                        false
                    } else {
                        true
                    }
                })
            }
        });

        log::info!("removeAlignables: {}", removed);
    }

    fn sort_fusion_results(&mut self) {
        self.m_fusion_results
            .sort_by(|a, b| Self::more_reads(a, b).unwrap())
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
}
