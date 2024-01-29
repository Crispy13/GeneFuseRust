use std::fmt::Write;
use std::io::Write as iowrite;
use std::{error::Error, io::BufWriter};

use crate::utils::StringCPP;

use super::sequence::{reverse_complement, Sequence};

#[derive(Debug, PartialEq, PartialOrd, Clone)]
pub(crate) struct SequenceRead {
    pub(crate) m_name: String,
    pub(crate) m_seq: Sequence,
    pub(crate) m_strand: String,
    pub(crate) m_quality: String,
    pub(crate) m_has_quality: bool,
}

impl SequenceRead {
    pub(crate) fn new(
        m_name: String,
        m_seq: String,
        m_strand: String,
        m_quality: String,
        m_has_quality: bool,
    ) -> Self {
        Self {
            m_name,
            m_seq: Sequence::new(m_seq),
            m_strand,
            m_quality,
            m_has_quality,
        }
    }

    pub(crate) fn len(&self) -> usize {
        self.m_seq.len()
    }

    fn print_with_breaks(&self, breaks: &[i32]) {
        // did not code this function at now.

        // std::cout << mName << endl;
        // std::cout << makeStringWithBreaks(mSeq.mStr, breaks)<< endl;
        // std::cout << mStrand << endl;
        // if(mHasQuality)
        //     std::cout << makeStringWithBreaks(mQuality, breaks) << endl;
    }

    pub(crate) fn print_html_td_with_breaks(
        &self,
        f: &mut std::io::BufWriter<std::fs::File>,
        breaks: Vec<i32>,
    ) -> Result<(), Box<dyn Error>> {
        write!(
            f,
            "<td class='alignright'>{}</td>",
            self.make_html_seq_with_qual(0, *breaks.get(0).unwrap())
        )?;

        for i in (0..(breaks.len() - 1)) {
            write!(f, "<td")?;
            if i == 0 {
                write!(f, " class='alignright'")?;
            }
            write!(
                f,
                ">{}</td>",
                self.make_html_seq_with_qual(
                    *breaks.get(i).unwrap(),
                    breaks.get(i + 1).unwrap() - breaks.get(i).unwrap()
                )
            )?;
        }

        if breaks.get(breaks.len() - 1).unwrap() > &0 {
            write!(
                f,
                "<td class='alignleft'>{}</td>",
                self.make_html_seq_with_qual(
                    *breaks.get(breaks.len() - 1).unwrap(),
                    self.m_seq.m_str.len() as i32 - *breaks.get(breaks.len() - 1).unwrap()
                )
            )?;
        }

        Ok(())
    }

    fn make_string_with_breaks(origin: &str, breaks: &[i32]) -> String {
        let mut ret = origin.subchars(0, *breaks.get(0).unwrap() as usize).to_string();

        for i in (0..(breaks.len() - 1)) {
            write!(
                &mut ret,
                " {}",
                origin.subchars(
                    *breaks.get(i).unwrap() as usize,
                    (breaks.get(i + 1).unwrap() - breaks.get(i).unwrap()) as usize
                )
            )
            .unwrap();
        }

        if *breaks.get(breaks.len() - 1).unwrap() > 0 {
            write!(
                &mut ret,
                " {}",
                origin.subchars(
                    *breaks.get(breaks.len() - 1).unwrap() as usize,
                    (origin.len() as i32 - breaks.get(breaks.len() - 1).unwrap()) as usize
                )
            )
            .unwrap();
        }

        ret
    }

    fn make_html_seq_with_qual(&self, start: i32, length: i32) -> String {
        let mut ss = String::new();
        for i in (start..(start + length).min(self.m_seq.len() as i32)).map(|e| e as usize) {
            write!(
                &mut ss,
                "<a title='{}'><font color='{}'>{}</font></a>",
                self.m_quality.get(i..(i+1)).unwrap(),
                quality_color(*self.m_quality.as_bytes().get(i).unwrap() as char),
                self.m_seq.m_str.get(i..(i+1)).unwrap(),
            )
            .unwrap();
        }

        ss
    }

    pub(crate) fn last_index(&self) -> String {
        let len = self.m_name.len();

        if len < 5 {
            return "".to_string();
        }

        for (ch, i) in self.m_name.chars().rev().zip((0..=(len - 1)).rev()).skip(4) {
            if ch == ':' || ch == '+' {
                return self.m_name.subchars(i + 1, len - i).to_string();
            }
        }

        "".to_string()
    }

    /// qual = 20 default.
    fn low_qual_count(&self, qual: i32) -> i32 {
        let mut count = 0;
        for q in self.m_quality.chars() {
            if (q as u32 as i32) < qual + 33 {
                count += 1;
            }
        }

        count
    }

    pub(crate) fn reverse_complement(&self) -> SequenceRead {
        let seq = reverse_complement(&self.m_seq.m_str);
        let qual = {
            let mut qual = self.m_quality.clone().into_bytes();
            qual.reverse();

            String::from_utf8(qual).unwrap()
        };
        let strand = {
            if self.m_strand == "+" {
                "-"
            } else {
                "+"
            }
        }
        .to_string();

        SequenceRead::new(self.m_name.clone(), seq, strand, qual, true)
    }

    pub(crate) fn print_file(
        &self,
        f: &mut BufWriter<std::fs::File>,
    ) -> Result<(), Box<dyn Error>> {
        writeln!(f, "{}", self.m_name)?;
        writeln!(f, "{}", self.m_seq.m_str)?;
        writeln!(f, "{}", self.m_strand)?;
        if self.m_has_quality {
            writeln!(f, "{}", self.m_quality)?;
        }

        Ok(())
    }
}

fn quality_color(qual: char) -> &'static str {
    if qual >= 'I' {
        // >= Q40, extremely high quality
        return "#78C6B9";
    }
    if qual >= '?' {
        // Q30 ~ Q39, high quality
        return "#33BBE2";
    }

    if qual >= '5' {
        // Q20 ~ Q29, moderate quality
        return "#666666";
    }

    if qual >= '0' {
        // Q15 ~ Q19, low quality
        return "#E99E5B";
    } else {
        // <= Q14, extremely low quality
        return "#FF0000";
    }
}

#[derive(Clone, Debug)]
pub(crate) struct SequenceReadPair {
    pub(crate) m_left: SequenceRead,
    pub(crate) m_right: SequenceRead,
}

impl SequenceReadPair {
    pub(crate) fn new(left: SequenceRead, right: SequenceRead) -> SequenceReadPair {
        Self {
            m_left: left,
            m_right: right,
        }
    }

    pub(crate) fn fast_merge(&self) -> Option<SequenceRead> {
        let rc_right = self.m_right.reverse_complement();
        let len1 = self.m_left.len();
        let len2 = self.m_right.len();

        // use the pointer directly for speed
        let str1 = self.m_left.m_seq.m_str.as_str();
        let str2 = rc_right.m_seq.m_str.as_str();
        let qual1 = self.m_left.m_quality.as_str();
        let qual2 = rc_right.m_quality.as_str();

        // we require at least 30 bp overlapping to merge a pair
        const MIN_OVERLAP: i32 = 30;

        let mut overlapped = false;
        let mut olen = MIN_OVERLAP;
        let mut diff = 0;

        // the diff count for 1 high qual + 1 low qual
        let mut low_qual_diff = 0;

        let str1_cv = str1.chars().collect::<Vec<char>>();
        let str2_cv = str2.chars().collect::<Vec<char>>();
        let qual1_cv = qual1.chars().collect::<Vec<char>>();
        let qual2_cv = qual2.chars().collect::<Vec<char>>();

        while olen <= len1.min(len2) as i32 {
            diff = 0;
            low_qual_diff = 0;
            let mut ok = true;
            let offset = len1 as i32 - olen;
            for i in (0..(olen)) {
                if str1_cv.get((offset + i) as usize).unwrap() != str2_cv.get(i as usize).unwrap() {
                    diff += 1;
                    // one is >= Q30 and the other is <= Q15
                    if (qual1_cv.get((offset + i) as usize).unwrap() >= &'?'
                        && qual2_cv.get(i as usize).unwrap() <= &'0')
                        || (qual1_cv.get((offset + i) as usize).unwrap() <= &'0'
                            && qual2_cv.get(i as usize).unwrap() >= &'?')
                    {
                        low_qual_diff += 1;
                    }
                    // we disallow high quality diff, and only allow up to 3 low qual diff
                    if diff > low_qual_diff || low_qual_diff >= 3 {
                        ok = false;
                        break;
                    }
                }
            }
            if ok {
                overlapped = true;
                break;
            }
            olen += 1;
        }

        if overlapped {
            let offset = len1 as i32 - olen;
            let mut ss = String::new();
            write!(&mut ss, "{} merged_diff_{}", self.m_left.m_name, diff).unwrap();
            let merged_name = ss;
            // let merged_seq = format!(
            //     "{}{}",
            //     self.m_left.m_seq.m_str.subchars(0, offset as usize),
            //     rc_right.m_seq.m_str
            // );
            let mut merged_seq_cv = str1_cv
                .get(..(offset as usize))
                .unwrap()
                .iter()
                .copied()
                .chain(str2_cv.iter().copied())
                .collect::<Vec<char>>();

            // let merged_qual = format!(
            //     "{}{}",
            //     self.m_left.m_quality.subchars(0, offset as usize),
            //     rc_right.m_quality
            // );
            let mut merged_qual_cv = qual1_cv
                .get(..(offset as usize))
                .unwrap()
                .iter()
                .copied()
                .chain(qual2_cv.iter().copied())
                .collect::<Vec<char>>();

            // quality adjuction and correction for low qual diff
            for i in (0..(olen)) {
                if str1_cv.get((offset + i) as usize).unwrap() != str2_cv.get(i as usize).unwrap() {
                    if qual1_cv.get((offset + i) as usize).unwrap() >= &'?'
                        && qual2_cv.get(i as usize).unwrap() <= &'0'
                    {
                        *merged_seq_cv.get_mut((offset + i) as usize).unwrap() =
                            *str1_cv.get((offset + i) as usize).unwrap();
                        *merged_qual_cv.get_mut((offset + i) as usize).unwrap() =
                            *qual1_cv.get((offset + i) as usize).unwrap()
                    } else {
                        *merged_seq_cv.get_mut((offset + i) as usize).unwrap() =
                            *str2_cv.get((i) as usize).unwrap();
                        *merged_qual_cv.get_mut((offset + i) as usize).unwrap() =
                            *qual2_cv.get((i) as usize).unwrap()
                    }
                } else {
                    // add the quality of the pair to make a high qual
                    *merged_qual_cv.get_mut((offset + i) as usize).unwrap() = char::from_u32(
                        *qual1_cv.get((offset + i) as usize).unwrap() as u32
                            + *qual2_cv.get(i as usize).unwrap() as u32
                            - 33,
                    )
                    .unwrap();
                    if merged_qual_cv.get((offset + i) as usize).unwrap() >= &'Z' {
                        *merged_qual_cv.get_mut((offset + i) as usize).unwrap() = 'Z'
                    }
                }
            }

            return Some(SequenceRead::new(
                merged_name,
                merged_seq_cv.into_iter().collect::<String>(),
                "+".to_string(),
                merged_qual_cv.into_iter().collect::<String>(),
                true,
            ));
        }

        None
    }
}

#[cfg(test)]
mod test {

    use crate::core::read::SequenceReadPair;

    use super::SequenceRead;

    fn _fast_merge() -> bool {
        let left = SequenceRead::new(
            "@NS500713:64:HFKJJBGXY:1:11101:20469:1097 1:N:0:TATAGCCT+GGTCCCGA".to_string(),
            "TTTTTTCTCTTGGACTCTAACACTGTTTTTTCTTATGAAAACACAGGAGTGATGACTAGTTGAGTGCATTCTTATGAGACTCATAGTCATTCTATGATGTAG".to_string(),
            "+".to_string(),
            "AAAAA6EEEEEEEEEEEEEEEEE#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE".to_string(),
            true
        );

        let right = SequenceRead::new(
            "@NS500713:64:HFKJJBGXY:1:11101:20469:1097 1:N:0:TATAGCCT+GGTCCCGA".to_string(),
            "AAAAAACTACACCATAGAATGACTATGAGTCTCATAAGAATGCACTCAACTAGTCATCACTCCTGTGTTTTCATAAGAAAAAACAGTGTTAGAGTCCAAGAG".to_string(),
            "+".to_string(),
            "AAAAA6EEEEE/EEEEEEEEEEE#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE".to_string(),
            true
        );

        let pair = SequenceReadPair::new(left, right);

        let merged = pair.fast_merge();

        if merged.is_none() {
            return false;
        }

        if merged.as_ref().unwrap().m_seq.m_str != "TTTTTTCTCTTGGACTCTAACACTGTTTTTTCTTATGAAAACACAGGAGTGATGACTAGTTGAGTGCATTCTTATGAGACTCATAGTCATTCTATGATGTAGTTTTTT"{
            return false;
        }

        println!("{:?}", merged);
        true
    }

    #[test]
    fn fast_merge() {
        assert_eq!(true, _fast_merge());
    }
}
