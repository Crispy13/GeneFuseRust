use crate::{
    core::{
        edit_distance::{self, edit_distance, edit_distance_from_str},
        overlap,
    },
    utils::StringCPP,
};

use super::sequence::Sequence;

pub(crate) struct Overlap {
    pub(crate) m_offset: i32,
    pub(crate) m_overlap_len: i32,
    pub(crate) m_distance: i32,
    pub(crate) m_overlapped: bool,
}

impl Overlap {
    pub(crate) fn new(offset: i32, overlap_len: i32, distance: i32) -> Self {
        Self {
            m_offset: offset,
            m_overlap_len: overlap_len,
            m_distance: distance,
            m_overlapped: overlap_len > 0,
        }
    }

    pub(crate) fn fit(r1: &Sequence, r2: &Sequence) -> Self {
        let len1 = r1.len() as i32;
        let len2 = r2.len() as i32;

        let reverse_r2 = r2.reverse_complement();

        let mut overlapped = false;
        let mut overlap_len = 0;
        let mut offset = 0_i32;
        let mut distance = 0_i32;

        // a match of less than 10 is considered as unconfident
        while (offset < len1 - 10 && overlapped == false) {
            // the overlap length of r1 & r2 when r2 is move right for offset
            overlap_len = (len1 - offset).min(len2);

            log::debug!("offset={}, overlap_len={}, r1={}, reversed_r2={}", offset, overlap_len, r1.m_str, reverse_r2.m_str);
            distance = edit_distance_from_str(
                &r1.m_str.subchars(offset as usize, overlap_len as usize),
                &reverse_r2.m_str.subchars(0, overlap_len as usize),
            ) as i32;
            log::debug!("distance={}", distance);

            let threshold = (3.0_f32).min(overlap_len as f32 / 10.0);
            log::debug!("threshold={}", threshold);
            if distance as f32 <= threshold {
                // now we find a good candidate
                // we verify it by moving r2 one more base to see if the distance is getting longer
                // if yes, then current is the best match, otherwise it's not
                while offset < len1 - 10 {
                    let next_offset = offset + 1;
                    let next_overlap_len = (len1 - next_offset).min(len2);
                    let next_distance = edit_distance_from_str(
                        &r1.m_str
                            .subchars(next_offset as usize, next_overlap_len as usize),
                        &reverse_r2.m_str.subchars(0, next_overlap_len as usize),
                    ) as i32;

                    log::debug!("next_offset={}, next_overlap_len={}", next_offset, next_overlap_len, );
                    log::debug!("next_distance={}", next_distance);

                    if distance <= next_distance {
                        overlapped = true;
                        break;
                    } else {
                        offset = next_offset;
                        distance = next_distance;
                        overlap_len = next_overlap_len;
                    }
                }
                break;
            } else {
                offset += 1.max((distance - (threshold.ceil() as i32)) / 2)
            }
        }

        if overlapped && offset == 0 {
            // check if distance can get smaller if offset goes negative
            // this only happens when insert DNA is shorter than sequencing read length, and some adapter/primer is sequenced but not trimmed cleanly
            // we go reversely
            while offset > -(len2 - 10) {
                // the overlap length of r1 & r2 when r2 is move right for offset
                overlap_len = len1.min(len2 - offset.abs());
                distance = edit_distance_from_str(
                    &r1.m_str.subchars(0, overlap_len as usize),
                    &reverse_r2
                        .m_str
                        .subchars(-offset as usize, overlap_len as usize),
                ) as i32;
                let threshold = 3.0_f32.min(overlap_len as f32 / 10.0);
                if distance as f32 <= threshold {
                    while offset > -(len2 - 10) {
                        let next_offset = offset - 1;
                        let next_overlap_len = len1.min(len2 - next_offset.abs());
                        let next_distance = edit_distance_from_str(
                            &r1.m_str.subchars(0, next_overlap_len as usize),
                            &reverse_r2
                                .m_str
                                .subchars(-next_offset as usize, next_overlap_len as usize),
                        ) as i32;
                        if distance <= next_distance {
                            return Self::new(offset, overlap_len, distance);
                        } else {
                            distance = next_distance;
                            overlap_len = next_overlap_len;
                            offset = next_offset;
                        }
                    }
                } else {
                    offset -= 1.max(distance - threshold.ceil() as i32 / 2);
                }
            }
        } else if overlapped {
            return Self::new(offset, overlap_len, distance);
        }

        Self::new(0, 0, 0)
    }
}

#[cfg(test)]
mod test {
    use crate::{core::{overlap::Overlap, sequence::Sequence}, utils::logging::init_logger};

    use std::sync::Once;

    

    #[test]
    fn minus() {
        let a = 2;
        println!("{}", -a as i32);
        println!("{}", (-a) as i32);
    }

    fn _overlap_test() -> bool {
        let r1 = [
            Sequence::new("TTTGCAGGCACCTACCACTGTACCTGTCTAATTTTTCTTCTGCCCTTTTTTTTTTTTTTTTTTTTTTTTTGGGGTAGAGACGAGGCCTTGCTATGTAGCCCTTGCTGGTCTCAAACTCCTCGCCTCAAGTGATCCTCCTGCCTCGGCCTCC".to_string()),
                Sequence::new("CCCTATGTCTACAAAACATCAGAAAATTAGGGTGTGGTGGCTCATGCCTATAGTCATAGCTACATAGGAGGCTGAGGCAGGAGGATCGCTTGAGGGCAGGAGGATCACTCGAGCTCTGAAGGTCAACGCTGCAGTGAGCTATGATCGTGCC".to_string()),
            Sequence::new("TAGAGGGCTCAGATGCATTCCTTTTTAGCAGTGCTCTTATTTGGCATTGGTGGTGCTGTTTCTGTTGACCACTCCCAGAGTCTCTGGATGTTTTGTTATTCCTTTACCTCCCTAGCCTCTCCTTGGGGTTTCTTTGCAGGCTCTTGCTCTC".to_string()),
            Sequence::new("CCTGGGTAGCTGGGATACAGGCGCCCGCCACCACGCCCGGCTAATTTTGTATTTTTAGTAGAGACGAGGTTTCACCACATTGGCCAGGCTGGTCTCAAACTCCTGACCTCAGGTGATCTGCCTGCCTCAGCCTCCTAGAGTGCTGGG".to_string()),
            Sequence::new("GTTCCTTTTAACATAGAAAGCAGCTAATTTTCCTATTCAAAAAATGGAGCTCTATTAAAAGATAAAACAGCAGCTTAGCTCTAGGTAAAGTGATCCATGCGGTTCTTCTTCTTTTTTTTGTTTTGAGATGGACTCTCGCTCTGTCACCCA".to_string()),
            
        ];

        let r2 = [
             Sequence::new("CATGGTGGCTCATGCCTGTAATCCCAGTGGTTTGGGAGGCCGAGGCAGGAGGATCACTTGAGGCGAGGAGTTTGAGACCAGCAAGGGCTACATAGCAAGGCCTCGTCTCTACCCCAAAAAAAAAAAAAAAAAAAAAAAAAGGGCAGAAGAA".to_string()),
            Sequence::new("AGTGCAGTGGCACGATCATAGCTCACTGCAGCGTTGACCTTCAGAGCTCGAGTGATCCTCCTGCCCTCAAGCGATCCTCCTGCCTCAGCCTCCTATGTAGCTATGACTATAGGCATGAGCCACCACACCCTAATTTTCTGATGTTTTGTAG".to_string()),
            Sequence::new("CTGGAGATAAACACCTAGCAGTCATGAGACAAAGCTCTGCAATGCTTGTATTTATGGGATACAAGAGAGAGCAAGAGCCTGCAAAGAAACCCCAAGGAGAGGCTAGGGAGGTAAAGGAATAACAAAACATCCAGAGACACTGGGAGTGGTC".to_string()),
            Sequence::new("CCCAGCACTCTAGGAGGCTGAGGCAGGCAGATCACCTGAGGTCAGGAGTTTGAGACCAGCCTGGCCAATGTGGTGAAACCTCGTCTCTACTAAAAATACAAAATTAGCCGGGCGTGGTGGCGGGCGCCTGTAATCCCAGCTACCCAGC".to_string()),
            Sequence::new("TGGGTGACAGAGCGAGAGTCCATCTCAAAACAAAAAAAAGAAGAAGAACCGCACTGGATCACTTTACCTCAGAGCTAAGCTGCTGTTTTATCTTTTAATAGAGCTCCATTTTTTGAATAGGAAAATTAGCTGCTTTCTATGTTAAAAGGAA".to_string())
        ];

        let overlap = [
            Overlap::new(34, 117, 0),
            Overlap::new(8, 143, 0),
            Overlap::new(66, 85, 1),
            Overlap::new(-1, 147, 2),
            Overlap::new(0, 0, 0)
        ];
        
        for i in (0..5) {
            let fit = Overlap::fit(&r1[i], &r2[i]);
            if fit.m_offset != overlap[i].m_offset || fit.m_overlap_len!=overlap[i].m_overlap_len || fit.m_distance!=overlap[i].m_distance {
                eprintln!("Fail in Overlap::fit() with sequence");
                eprintln!("{}", r1[i].m_str);
                eprintln!("{}", r2[i].m_str);
                return false;
            }
        }

        true


    }

    #[test]
    fn overlap_test() {
        // set_logger();

        assert_eq!(true, _overlap_test());
    }

    #[test]
    fn cast_f64() {
        println!("{}", (3.0_f32).min(151 as f32 / 10.0));
        println!("{}", 1.max(68 - (3_f32.ceil() as i32) / 2));
    }
}
