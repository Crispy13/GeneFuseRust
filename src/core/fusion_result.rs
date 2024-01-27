use crate::{
    aux::global_settings::global_settings,
    core::{
        edit_distance::{edit_distance, edit_distance_from_str},
        sequence::{reverse_complement, Sequence},
    },
    utils::{dis_connected_count, int2str, StringCPP},
};

use super::{
    common::GenePos, edit_distance, fusion::Fusion, gene::Gene, read, read_match::ReadMatch,
};
use std::{
    error::Error,
    fmt::Write,
    fs::File,
    io::{BufWriter, Read, Write as io_write},
    mem,
    ops::Neg,
};

#[derive(Default)]
pub(crate) struct FusionResult {
    pub(crate) m_left_gp: GenePos,
    pub(crate) m_right_gp: GenePos,
    pub(crate) m_matches: Vec<ReadMatch>,
    pub(crate) m_unique: i32,
    pub(crate) m_title: String,
    pub(crate) m_left_ref: String,
    pub(crate) m_right_ref: String,
    pub(crate) m_left_ref_ext: String,
    pub(crate) m_right_ref_ext: String,
    pub(crate) m_left_pos: String,
    pub(crate) m_right_pos: String,
    pub(crate) m_left_gene: Gene,
    pub(crate) m_right_gene: Gene,
    pub(crate) m_left_is_exon: bool,
    pub(crate) m_right_is_exon: bool,
    pub(crate) m_left_exon_or_intron_id: i32,
    pub(crate) m_right_exon_or_intron_id: i32,
    pub(crate) m_left_exon_num: f32,
    pub(crate) m_left_intron_num: f32,
    pub(crate) m_right_exon_num: f32,
    pub(crate) m_right_intron_num: f32,
}

impl FusionResult {
    pub(crate) fn with_minimum() -> Self {
        Self {
            m_left_is_exon: false,
            m_right_is_exon: false,
            m_left_exon_or_intron_id: -1,
            m_right_exon_or_intron_id: -1,
            ..Default::default()
        }
    }

    pub(crate) fn calc_fusion_point(&mut self) -> () {
        if self.m_matches.len() == 0 {
            return;
        }

        // if we can find an exact match with 0 gap, then use it
        // else we use the mean one
        let mut left_total = 0;
        let mut right_total = 0;

        for read_match in self.m_matches.iter() {
            if read_match.m_gap == 0 {
                self.m_left_gp = read_match.m_left_gp.clone();
                self.m_right_gp = read_match.m_right_gp.clone();

                return;
            }

            left_total += read_match.m_left_gp.position;
            right_total += read_match.m_right_gp.position;
        }

        self.m_left_gp.contig = self.m_matches.first().unwrap().m_left_gp.contig;
        self.m_left_gp.contig = (left_total / self.m_matches.len() as i32) as i16;
        self.m_right_gp.contig = self.m_matches.first().unwrap().m_right_gp.contig;
        self.m_right_gp.contig = (right_total / self.m_matches.len() as i32) as i16;
    }

    pub(crate) fn calc_unique(&mut self) -> () {
        self.m_unique = 1;

        if self.m_matches.len() < 2 {
            return;
        }

        // since it is sorted, so just check every match with previous one
        let mut mm_iter = self.m_matches.iter();
        let mut prev_mm = mm_iter.next().unwrap();

        while let Some(mm) = mm_iter.next() {
            if mm.m_read_break != prev_mm.m_read_break || mm.m_read.len() != prev_mm.m_read.len() {
                self.m_unique += 1;
            }
            prev_mm = mm;
        }
    }

    pub(crate) fn is_deletion(&self) -> bool {
        if self.m_left_gp.contig == self.m_right_gp.contig {
            if self.m_left_gp.position > 0 && self.m_right_gp.position > 0 {
                return true;
            }
            if self.m_left_gp.position < 0 && self.m_right_gp.position < 0 {
                return true;
            }
        }

        false
    }

    fn can_be_mapped(&self) -> bool {
        if self.can_be_matched(&self.m_left_ref_ext, &self.m_right_ref) {
            return true;
        }
        if self.can_be_matched(&self.m_left_ref, &self.m_right_ref_ext) {
            return true;
        }

        false
    }

    fn can_be_matched(&self, s1: &str, s2: &str) -> bool {
        let len = s1.len();
        for offset in (-6..=6_i32) {
            let start1 = offset.max(0);
            let start2 = offset.neg().max(0);
            let cmplen = len as i32 - offset.abs();
            if start1 >= (s1.len() as i32) || start2 >= (s2.len() as i32) {
                return true;
            }
            let ed = edit_distance_from_str(
                &s1.subchars(start1 as usize, cmplen as usize),
                &s2.subchars(start2 as usize, cmplen as usize),
            );

            let threshold = cmplen / 10;
            if (ed as i32) <= threshold {
                return true;
            }
        }

        false
    }

    pub(crate) fn is_qualified(&self) -> bool {
        if self.m_unique < global_settings().unique_requirement as i32 {
            return false;
        }

        if self.can_be_mapped() {
            return false;
        }

        if self.m_left_ref.len() <= 30 || self.m_right_ref.len() <= 30 {
            return false;
        }

        if dis_connected_count(&self.m_left_ref.subchars(self.m_left_ref.len() - 10, 10)) <= 2 {
            return false;
        }

        if dis_connected_count(&self.m_right_ref.subchars(0, 10)) <= 2 {
            return false;
        }

        true
    }

    pub(crate) fn update_info(&mut self, fusions: &[Fusion]) -> () {
        self.m_left_gene = fusions
            .get(self.m_left_gp.contig as usize)
            .unwrap()
            .m_gene
            .clone();
        self.m_right_gene = fusions
            .get(self.m_right_gp.contig as usize)
            .unwrap()
            .m_gene
            .clone();

        let mut ss = String::new();
        if self.is_deletion() {
            ss.push_str("Deletion: ");
        } else {
            ss.push_str("Fusion: ");
        }

        write!(
            ss,
            "{}___{}  (total: {}, unique:{})",
            self.m_left_gene.pos2str(self.m_left_gp.position).unwrap(),
            self.m_right_gene.pos2str(self.m_right_gp.position).unwrap(),
            self.m_matches.len(),
            self.m_unique
        )
        .unwrap();

        self.m_title = ss;

        self.m_left_pos = self.m_left_gene.pos2str(self.m_left_gp.position).unwrap();
        self.m_right_pos = self.m_right_gene.pos2str(self.m_right_gp.position).unwrap();

        self.m_left_gene.get_exon_intron(
            self.m_left_gp.position,
            self.m_left_is_exon,
            self.m_left_exon_or_intron_id,
        );
        self.m_right_gene.get_exon_intron(
            self.m_right_gp.position,
            self.m_right_is_exon,
            self.m_right_exon_or_intron_id,
        );
    }

    pub(crate) fn make_reference(&mut self, ref_l: &str, ref_r: &str) {
        let mut longest_left = 0;
        let mut longest_right = 0;

        for read_match in self.m_matches.iter() {
            if read_match.m_read_break + 1 > longest_left {
                longest_left = read_match.m_read_break + 1;
            }
            if read_match.m_read.len() - (read_match.m_read_break as usize + 1)
                > longest_right as usize
            {
                longest_right =
                    (read_match.m_read.len() - (read_match.m_read_break as usize + 1)) as i32;
            }
        }

        self.m_left_ref = get_ref_seq(
            ref_l,
            self.m_left_gp.position - longest_left + 1,
            self.m_left_gp.position,
        );
        self.m_right_ref = get_ref_seq(
            ref_r,
            self.m_right_gp.position,
            self.m_right_gp.position + longest_right - 1,
        );

        self.m_left_ref_ext = get_ref_seq(
            ref_l,
            self.m_left_gp.position,
            self.m_left_gp.position + longest_right - 1,
        );
        self.m_right_ref_ext = get_ref_seq(
            ref_r,
            self.m_right_gp.position - longest_left + 1,
            self.m_right_gp.position,
        );
    }

    pub(crate) fn adjust_fusion_break(&mut self) {
        let mut m_matches = mem::take(&mut self.m_matches);

        for (i, read_match) in m_matches.iter_mut().enumerate() {
            let mut smallest_ed = 0xFFFF;
            let mut shift = 0;

            for s in (-3..=3) {
                let mut left_ed = 0;
                let mut right_ed = 0;
                let ed = self.calc_ed(&read_match, s, &mut left_ed, &mut right_ed);
                if ed < smallest_ed {
                    smallest_ed = ed;
                    shift = s;
                    read_match.m_left_distance = left_ed;
                    read_match.m_right_distance = right_ed;
                }
            }

            read_match.m_read_break += shift;
            read_match.m_left_gp.position += shift;
            read_match.m_right_gp.position += shift;
        }

        self.m_matches = m_matches;
    }

    fn calc_ed(&self, m: &ReadMatch, shift: i32, left_ed: &mut i32, right_ed: &mut i32) -> i32 {
        let read_break = m.m_read_break + shift;
        let seq = m.m_read.m_seq.m_str.as_str();

        let left_len = read_break + 1;
        let right_len = seq.chars().count() as i32 - left_len;
        let left_seq = seq
            .chars()
            .skip(0)
            .take(left_len as usize)
            .collect::<String>();
        let right_seq = seq
            .chars()
            .skip(left_len as usize)
            .take(right_len as usize)
            .collect::<String>();

        // use the sequence near the break point to adjust
        let mut left_comp = left_seq.len().min(self.m_left_ref.len()).min(20) as i32;
        let mut right_comp = right_seq.len().min(self.m_right_ref.len()).min(20) as i32;

        let left_part_ed = edit_distance_from_str(
            &left_seq
                .chars()
                .skip(left_seq.chars().count() - left_comp as usize)
                .take(left_comp as usize)
                .collect::<String>(),
            &self
                .m_left_ref
                .chars()
                .skip(self.m_left_ref.chars().count() - left_comp as usize)
                .take(left_comp as usize)
                .collect::<String>(),
        );

        let right_part_ed = edit_distance_from_str(
            &right_seq
                .chars()
                .skip(0)
                .take(right_comp as usize)
                .collect::<String>(),
            &self
                .m_right_ref
                .chars()
                .skip(0)
                .take(right_comp as usize)
                .collect::<String>(),
        );

        let total_ed = left_part_ed + right_part_ed;

        // recalculate the left and right edit distance, but
        left_comp = left_len.min(self.m_left_ref.chars().count() as i32) as i32;
        right_comp = right_len.min(self.m_right_ref.chars().count() as i32) as i32;

        *left_ed = edit_distance_from_str(
            &left_seq
                .chars()
                .skip(left_seq.chars().count() - left_comp as usize)
                .take(left_comp as usize)
                .collect::<String>(),
            &self
                .m_left_ref
                .chars()
                .skip(self.m_left_ref.chars().count() - left_comp as usize)
                .take(left_comp as usize)
                .collect::<String>(),
        ) as i32;

        *right_ed = edit_distance_from_str(
            &right_seq
                .chars()
                .skip(0)
                .take(right_comp as usize)
                .collect::<String>(),
            &self
                .m_right_ref
                .chars()
                .skip(0)
                .take(right_comp as usize)
                .collect::<String>(),
        ) as i32;

        total_ed as i32
    }

    pub(crate) fn add_match(&mut self, m: ReadMatch) {
        self.m_matches.push(m);
    }

    pub(crate) fn support(&self, m: &ReadMatch) -> bool {
        for read_match in self.m_matches.iter() {
            if self.support_same(m, read_match) {
                return true;
            }
        }

        false
    }

    fn support_same(&self, m1: &ReadMatch, m2: &ReadMatch) -> bool {
        const T: i32 = 3;

        if (m1.m_left_gp.position - m2.m_left_gp.position).abs() > T {
            return false;
        }
        if (m1.m_right_gp.position - m2.m_right_gp.position).abs() > T {
            return false;
        }

        if (m1.m_left_gp.contig != m2.m_left_gp.contig) {
            return false;
        }

        if (m1.m_right_gp.contig != m2.m_right_gp.contig) {
            return false;
        }

        true
    }
    pub(crate) fn is_left_protein_forward(&self) -> bool {
        if self.m_left_gene.is_reversed() {
            self.m_left_gp.position < 0
        } else {
            self.m_left_gp.position > 0
        }
    }

    pub(crate) fn is_right_protein_forward(&self) -> bool {
        if self.m_right_gene.is_reversed() {
            self.m_right_gp.position < 0
        } else {
            self.m_right_gp.position > 0
        }
    }

    fn calc_left_exon_intron_number(&mut self) -> () {
        let total_exon = self.m_left_gene.m_exons.len();
        let total_intron = total_exon - 1;

        if self.is_left_protein_forward() {
            if self.m_left_is_exon {
                self.m_left_exon_num = self.m_left_exon_or_intron_id as f32 - 0.5;
                self.m_left_intron_num = self.m_left_exon_or_intron_id as f32 - 1.0;
            } else {
                self.m_left_exon_num = self.m_left_exon_or_intron_id as f32;
                self.m_left_intron_num = self.m_left_exon_or_intron_id as f32 - 0.5;
            }
        } else {
            if self.m_left_is_exon {
                self.m_left_exon_num =
                    (total_exon as i32 - self.m_left_exon_or_intron_id) as f32 + 0.5;
                self.m_left_intron_num =
                    (total_intron as i32 - self.m_left_exon_or_intron_id) as f32 + 1.0;
            } else {
                self.m_left_exon_num = (total_exon as i32 - self.m_left_exon_or_intron_id) as f32;
                self.m_left_intron_num =
                    (total_intron as i32 - self.m_left_exon_or_intron_id) as f32 + 0.5;
            }
        }
    }

    fn calc_right_exon_intron_number(&mut self) -> () {
        let total_exon = self.m_right_gene.m_exons.len();
        let total_intron = total_exon - 1;

        if self.is_right_protein_forward() {
            if self.m_right_is_exon {
                self.m_right_exon_num =
                    (total_exon as i32 - self.m_right_exon_or_intron_id) as f32 + 0.5;
                self.m_right_intron_num =
                    (total_intron as i32 - self.m_right_exon_or_intron_id) as f32 + 1.0;
            } else {
                self.m_right_exon_num = (total_exon as i32 - self.m_right_exon_or_intron_id) as f32;
                self.m_right_intron_num =
                    (total_intron as i32 - self.m_right_exon_or_intron_id) as f32 + 0.5;
            }
        } else {
            if self.m_right_is_exon {
                self.m_right_exon_num = self.m_right_exon_or_intron_id as f32 - 0.5;
                self.m_right_intron_num = self.m_right_exon_or_intron_id as f32 - 1.0;
            } else {
                self.m_right_exon_num = self.m_right_exon_or_intron_id as f32;
                self.m_right_intron_num = self.m_right_exon_or_intron_id as f32 - 0.5;
            }
        }
    }

    pub(crate) fn print_fusion_protein_html(
        &mut self,
        buf_writer: &mut BufWriter<File>,
    ) -> Result<(), Box<dyn Error>> {
        self.calc_left_exon_intron_number();
        self.calc_right_exon_intron_number();

        let left_size = self.m_left_exon_num + self.m_left_intron_num;
        let right_size = self.m_right_exon_num + self.m_right_intron_num;

        let mut left_percent = (left_size * 100.0 / (left_size + right_size)).round() as i32;
        let mut right_percent = 100 - left_percent;
        if left_percent == 0 {
            left_percent = 1;
        }
        if right_percent == 0 {
            right_percent = 1;
        }

        write!(
            buf_writer,
            "{}",
            "<table width='100%' class='protein_table'>\n"
        )?;
        write!(buf_writer, "{}", "<tr>")?;
        write!(buf_writer, "<td width='{}%'>", left_percent)?;
        write!(buf_writer, "{}", self.m_left_gene.m_name)?;
        write!(buf_writer, "{}", "</td>")?;
        write!(buf_writer, "<td width='{}%'>", right_percent)?;
        write!(buf_writer, "{}", self.m_right_gene.m_name)?;
        write!(buf_writer, "{}", "</td>")?;
        write!(buf_writer, "{}", "</tr>")?;
        write!(buf_writer, "{}", "<tr>")?;
        write!(
            buf_writer,
            "<td class='protein_left' width='{}%'>",
            left_percent
        )?;

        self.print_left_protein_html(buf_writer)?;
        write!(buf_writer, "{}", "</td>")?;
        write!(
            buf_writer,
            "<td class='protein_right' width='{}%'>",
            left_percent
        )?;

        self.print_right_protein_html(buf_writer)?;
        write!(buf_writer, "{}", "</td>")?;
        write!(buf_writer, "{}", "</tr>")?;
        write!(buf_writer, "{}", "</table>")?;

        Ok(())
    }

    fn print_left_protein_html(
        &self,
        buf_writer: &mut BufWriter<File>,
    ) -> Result<(), Box<dyn Error>> {
        let total_step = self.m_left_exon_num + self.m_left_intron_num;
        let mut exon = 1_i32;
        let mut intron = 1;
        let mut step = 1_i32;

        let step_percent = 100.0 / total_step;
        let half_step_percent = step_percent * 0.5;
        let forward = self.is_left_protein_forward();

        if !forward {
            exon = self.m_left_gene.m_exons.len() as i32;
            intron = exon - 1;
            step = -1;
        }

        write!(
            buf_writer,
            "<table width='100%' class='protein_table'>\n<tr>"
        )?;

        let mut print_exon = 0_f32;
        let mut print_intron = 0_f32;

        while (print_exon < self.m_left_exon_num) || print_intron < self.m_left_intron_num {
            if print_exon < self.m_left_exon_num {
                let mut percent = step_percent;
                // last one is a half exon
                if print_exon + 1.0 > self.m_left_exon_num {
                    percent = half_step_percent;
                }
                self.print_exon_intron_td(
                    buf_writer,
                    true,
                    forward,
                    exon as i32,
                    percent,
                    "exon_left",
                )?;
                print_exon += 1.0;
                exon += step;
            }

            if print_intron < self.m_left_intron_num {
                let mut percent = step_percent;
                // last one is a half intron
                if print_intron + 1.0 > self.m_left_intron_num {
                    percent = half_step_percent;
                }
                self.print_exon_intron_td(
                    buf_writer,
                    false,
                    forward,
                    intron as i32,
                    percent,
                    "intron_left",
                )?;

                print_intron += 1.0;
                intron += step;
            }
        }

        write!(buf_writer, "</tr></table>")?;

        Ok(())
    }

    fn print_right_protein_html(
        &self,
        buf_writer: &mut BufWriter<File>,
    ) -> Result<(), Box<dyn Error>> {
        let total_step = self.m_right_exon_num + self.m_right_intron_num;
        let mut exon = self.m_right_exon_or_intron_id;
        let mut intron = self.m_right_exon_or_intron_id;
        let mut step = 1_i32;

        let step_percent = 100.0 / total_step;
        let half_step_percent = step_percent * 0.5;
        let forward = self.is_right_protein_forward();

        if !forward {
            step = -1;
        }

        write!(
            buf_writer,
            "<table width='100%' class='protein_table'>\n<tr>"
        )?;

        let mut print_exon = 0_f32;
        let mut print_intron = 0_f32;

        // print the first half intron
        if !self.m_right_is_exon {
            self.print_exon_intron_td(
                buf_writer,
                false,
                forward,
                intron,
                half_step_percent,
                "intron_right",
            )?;
            print_intron += 0.5;
            intron += step;
            if forward {
                exon += step;
            }
        }

        while (print_exon < self.m_right_exon_num) || print_intron < self.m_right_intron_num {
            if print_exon < self.m_right_exon_num {
                let mut percent = step_percent;
                if self.m_right_is_exon && print_exon == 0.0 {
                    percent = half_step_percent;
                }
                self.print_exon_intron_td(buf_writer, true, forward, exon, percent, "exon_right")?;
                if self.m_right_is_exon && print_exon == 0.0 {
                    print_exon += 0.5;
                } else {
                    print_exon += 1.0;
                }
                exon += step;
            }

            if print_intron < self.m_right_intron_num {
                let percent = step_percent;
                self.print_exon_intron_td(
                    buf_writer,
                    false,
                    forward,
                    intron,
                    percent,
                    "intron_right",
                )?;
            }
        }

        write!(buf_writer, "</tr></table>")?;

        Ok(())
    }

    fn print_exon_intron_td(
        &self,
        buf_writer: &mut BufWriter<File>,
        is_exon: bool,
        forward: bool,
        number: i32,
        percent: f32,
        style: &str,
    ) -> Result<(), Box<dyn Error>> {
        let mut int_percent = percent as i32;
        if int_percent <= 0 {
            int_percent = 1;
        }
        write!(
            buf_writer,
            "<td class='{}' width='{}%'>",
            style, int_percent
        )?;

        if is_exon {
            write!(buf_writer, "E{}", number)?;
        } else {
            if forward {
                write!(buf_writer, "→")?;
            } else {
                write!(buf_writer, "←")?;
            }
        }

        write!(buf_writer, "</td>")?;

        Ok(())
    }

    pub(crate) fn print(&self, fusion_list: &[Fusion]) {
        println!("\n#{}", self.m_title);
        for (i, m) in self.m_matches.iter().enumerate() {
            print!(">{}, ", i + 1);
            m.print();
        }
    }
}

fn get_ref_seq(ref_s: &str, start: i32, end: i32) -> String {
    // check start and end are in same strand
    if (start >= 0 && end <= 0) || (start <= 0 && end >= 0) {
        return "".into();
    }

    // check the overflow
    if start.abs() >= (ref_s.len() as i32) || end.abs() >= (ref_s.len() as i32) {
        return "".into();
    }

    let len = ((end - start).abs() + 1) as usize;

    if start < 0 {
        reverse_complement(
            &ref_s
                .chars()
                .skip(ref_s.chars().count() - end as usize)
                .take(len)
                .collect::<String>(),
        )
    } else {
        ref_s
            .chars()
            .skip(start as usize)
            .take(len)
            .collect::<String>()
    }
}

#[cfg(test)]
mod test {
    #[test]
    fn str_len() {
        let a = "abc";
        let b = "하이!";

        println!("{} {} {}", a.len(), b.len(), b.chars().count());

        // println!("{}", &b[..2]);
    }
}
