use core::fmt;
use std::{
    cmp::Ordering,
    error,
    fs::File,
    io::{BufWriter, Read, Write},
};

use crate::utils::StringCPP;

use super::{
    common::GenePos,
    fusion_scan::Error,
    read::{SequenceRead, SequenceReadCow, SequenceReadPair, SequenceReadPairCow},
};

#[derive(PartialEq, Clone, Debug)]
pub(crate) struct ReadMatch<'s> {
    pub(crate) m_read: SequenceRead,
    pub(crate) m_original_reads: Vec<SequenceReadCow<'s>>,
    pub(crate) m_overall_distance: i32,
    pub(crate) m_left_distance: i32,
    pub(crate) m_right_distance: i32,
    // the gap between left and right segment after segmentation
    pub(crate) m_gap: i32,
    pub(crate) m_reversed: bool,
    pub(crate) m_read_break: i32,
    pub(crate) m_left_gp: GenePos,
    pub(crate) m_right_gp: GenePos,
}

impl<'s> ReadMatch<'s> {
    /// reversed = false, on default
    pub(crate) fn new(
        r: SequenceRead,
        read_break: i32,
        left_gp: GenePos,
        right_gp: GenePos,
        gap: i32,
        reversed: bool,
    ) -> Self {
        Self {
            m_read: r,
            m_original_reads: Vec::new(),
            m_overall_distance: 0,
            m_left_distance: 0,
            m_right_distance: 0,
            m_gap: gap,
            m_reversed: reversed,
            m_read_break: read_break,
            m_left_gp: left_gp,
            m_right_gp: right_gp,
        }
    }

    pub(crate) fn less(m1: &ReadMatch, m2: &ReadMatch) -> bool {
        m1 < m2
    }

    pub(crate) fn greater(m1: &ReadMatch, m2: &ReadMatch) -> bool {
        m1 > m2
    }

    pub(crate) fn set_reversed(&mut self, flag: bool) {
        self.m_reversed = flag;
    }

    pub(crate) fn add_original_read(&mut self, r: SequenceReadCow<'s>) {
        self.m_original_reads.push(r);
    }

    pub(crate) fn add_original_pair(&mut self, pair: SequenceReadPairCow<'s>) -> () {
        let (p_m_left, p_m_right) = match pair {
            SequenceReadPairCow::Borrowed(b) => (
                SequenceReadCow::Borrowed(&b.m_left),
                SequenceReadCow::Borrowed(&b.m_right),
            ),
            SequenceReadPairCow::Owned(o) => (
                SequenceReadCow::Owned(o.m_left),
                SequenceReadCow::Owned(o.m_right),
            ),
        };

        self.m_original_reads.push(p_m_left);
        self.m_original_reads.push(p_m_right);
    }

    pub(crate) fn get_read(&self) -> &SequenceRead {
        &self.m_read
    }

    pub(crate) fn print_html_td(&self, f: &mut BufWriter<File>) -> Result<(), Error> {
        // write!(f, "d:{}", self.m_distance)?;

        if self.m_reversed {
            write!(f, "←")?;
        } else {
            write!(f, "→")?;
        }

        write!(f, "</a></span>")?;
        write!(
            f,
            "</td><td>{}|{}</td>",
            self.m_left_distance, self.m_right_distance
        )?;

        let mut breaks = Vec::new();
        breaks.push(self.m_read_break + 1);
        self.m_read.print_html_td_with_breaks(f, breaks)?;

        Ok(())
    }

    pub(crate) fn print_reads_to_file(&self, f: &mut BufWriter<File>) -> Result<(), Error> {
        self.m_original_reads
            .iter()
            .map(|r| r.print_file(f))
            .collect::<Result<(), Error>>()
    }

    pub(crate) fn print_read_to_json(
        &self,
        f: &mut BufWriter<File>,
        pad: &str,
    ) -> Result<(), Error> {
        writeln!(f, "{}\"seq\":\"{}\",", pad, self.m_read.m_seq.m_str)?;
        writeln!(f, "{}\"qual\":\"{}\"", pad, self.m_read.m_quality)?;

        Ok(())
    }

    pub(crate) fn print(&self) {
        print!("break:{}", self.m_read_break + 1);
        print!(
            ", diff:({} {})",
            self.m_left_distance, self.m_right_distance
        );

        if self.m_reversed {
            print!(", read direction: reversed complement");
        } else {
            print!(", read direction: original direction");
        }

        print!(
            ", name: {}",
            self.m_read.m_name.subchars(1, self.m_read.m_name.len() - 1)
        );
        print!("\n");
        print!(
            "{}",
            self.m_read
                .m_seq
                .m_str
                .subchars(0, (self.m_read_break + 1) as usize)
        );
        print!(" ");
        print!(
            "{}",
            self.m_read.m_seq.m_str.subchars(
                (self.m_read_break + 1) as usize,
                self.m_read.len() - (self.m_read_break + 1) as usize
            )
        );
        print!("\n");
    }
}

impl<'s> fmt::Display for ReadMatch<'s> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let dir_str = if self.m_reversed {
            ", read direction: reversed complement"
        } else {
            ", read direction: original direction"
        };

        write!(
            f,
            "break:{}, diff:({} {}){}, name: {}\n{} {}\n",
            self.m_read_break + 1,
            self.m_left_distance,
            self.m_right_distance,
            dir_str,
            self.m_read.m_name.chars().skip(1).collect::<String>(),
            self.m_read
                .m_seq
                .m_str
                .chars()
                .take(self.m_read_break as usize + 1)
                .collect::<String>(),
            self.m_read
                .m_seq
                .m_str
                .chars()
                .skip(self.m_read_break as usize + 1)
                .take(self.m_read.m_seq.m_str.chars().count() - (self.m_read_break as usize + 1))
                .collect::<String>(),
        )
    }
}

impl<'s> PartialOrd for ReadMatch<'s> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match self.m_read_break.partial_cmp(&other.m_read_break) {
            Some(Ordering::Equal) => {}
            ord => return ord,
        }

        // self.m_read
        //     .m_seq
        //     .m_str
        //     .len()
        //     .partial_cmp(&other.m_read.m_seq.m_str.len());

        match other
            .m_read
            .m_seq
            .m_str
            .len()
            .partial_cmp(&self.m_read.m_seq.m_str.len())
        {
            Some(Ordering::Equal) => {}
            ord => return ord,
        }

        self.m_read.m_name.partial_cmp(&other.m_read.m_name) // added 240201, not in original cpp code.
    }
}


#[cfg(test)]
mod test {
    use super::*;
    
    #[test]
    fn check_enums() {
        // let sp = SequenceReadPair::new(SequenceRead::new("ACGCGAGTA"), right)

        // let s = SequenceReadPairCow::Borrowed()
    }
}