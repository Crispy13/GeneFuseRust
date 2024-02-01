use std::collections::HashMap;
use std::hash::Hash;
use std::{fmt, iter};
use std::ops::BitAnd;

use super::matcher::GFHasherBuilder;

fn bp() {

}

fn edit_distance_bpv<T, const N: usize>(
    cmap: &mut HashMap<T, [u64; N], GFHasherBuilder>,
    v: &[T],
    vsize: usize,
    tmax: usize,
    tlen: usize,
) -> usize
where
    T: Hash,
    T: std::cmp::Eq,
    T: fmt::Debug,
    T: fmt::Display,
    T: Copy,
{
    let mut d = tmax * 64 + tlen;

    let top = 1_u64.wrapping_shl((tlen - 1) as u32);
    let lmb = 1_u64.wrapping_shl(63);

    // eprintln!("top={top}");
    // eprintln!("lmb={lmb}");
    // eprintln!("vsize={vsize}");

    let mut d0 = [0; N];
    let mut hp = [0; N];
    let mut hn = [0; N];

    let mut vp = vec![0_u64; tmax + 1];
    let mut vn = vec![0_u64; tmax + 1];

    (0..tmax).for_each(|i| *vp.get_mut(i).unwrap() = !0);
    (0..tlen).for_each(|i| {
        *vp.get_mut(tmax).unwrap() |= (1_u64.wrapping_shl(i as u32));
    });

    for i in (0..vsize) {
        let ch = v.get(i).unwrap();
        // let pm = match cmap.get_mut(ch) {
        //     Some(v) => v,
        //     None => panic!("cmap.keys()={:#?}, ch={}, i={}, v={:#?}",cmap.keys(), ch, i, v),
        // };
        let pm = cmap.entry(*ch).or_insert_with(|| [0; N]);

        for r in (0..=tmax) {
            let mut x = pm.get(r).unwrap().clone();
            if r > 0 && (hn.get(r - 1).unwrap() & lmb) != 0 {
                x |= 1_u64;
            }

            *d0.get_mut(r).unwrap() = (((x & vp.get(r).unwrap()).wrapping_add(*vp.get(r).unwrap()))
                ^ vp.get(r).unwrap())
                | x
                | vn.get(r).unwrap();

            *hp.get_mut(r).unwrap() =
                vn.get(r).unwrap() | !(d0.get(r).unwrap() | vp.get(r).unwrap());
            *hn.get_mut(r).unwrap() = d0.get(r).unwrap() & vp.get(r).unwrap();

            x = hp.get(r).unwrap().wrapping_shl(1);

            if (r == 0 || hp.get(r - 1).unwrap().bitand(lmb) != 0) {
                x |= 1_u64;
            }
            *vp.get_mut(r).unwrap() = (hn.get(r).unwrap().wrapping_shl(1)) | !(d0.get(r).unwrap() | x);
            if r > 0 && (hn.get(r - 1).unwrap() & lmb) != 0 {
                *vp.get_mut(r).unwrap() |= 1;
            }
            *vn.get_mut(r).unwrap() = d0.get(r).unwrap() & x;
        }

        // eprintln!("{}", hp.get(tmax).unwrap() & top);
        // eprintln!("{}", hn.get(tmax).unwrap() & top);
        if (hp.get(tmax).unwrap() & top) != 0 {
            d += 1;
        } else if (hn.get(tmax).unwrap() & top) != 0 {
            d -= 1;
        }
    }

    d
}

fn edit_distance_dp(str1: &str, size1: usize, str2: &str, size2: usize) -> usize {
    let mut d = iter::repeat_with(|| Vec::<u32>::with_capacity(size2 + 1))
        .take(size1 + 1)
        .collect::<Vec<_>>();

    (0..(size1 + 1)).for_each(|i| *d.get_mut(i).unwrap().get_mut(0).unwrap() = i as u32);
    (0..(size2 + 1)).for_each(|i| *d.get_mut(0).unwrap().get_mut(i).unwrap() = i as u32);

    for i in (1..(size1 + 1)) {
        let (d1p, d1) = d.split_at_mut(i);
        let d1p = d1p.last().unwrap();
        let d1 = d1.first_mut().unwrap();

        for j in (1..(size2 + 1)) {
            *d1.get_mut(j).unwrap() = (d1p.get(j).unwrap().min(d1.get(j - 1).unwrap()) + 1).min(
                d1p.get(j - 1).unwrap()
                    + if str1.chars().nth(i - 1).unwrap() == str2.chars().nth(j - 1).unwrap() {
                        0
                    } else {
                        1
                    },
            )
        }
    }

    d.get(size1).unwrap().get(size2).unwrap().clone() as usize
}

type VARR<const N: usize> = [u64; N];
fn edit_distance_map_<const N: usize>(a: &str, asize: usize, b: &str, bsize: usize) -> usize {
    let mut cmap: HashMap<char, VARR<N>, GFHasherBuilder> = HashMap::with_hasher(GFHasherBuilder::default());
    let tmax = (asize - 1).wrapping_shr(6);
    let tlen = asize - tmax * 64;

    // eprintln!("asize={:b}", asize-1);
    // eprintln!("tmax={}", tmax);
    // eprintln!("tlen={}", tlen);
    for i in (0..tmax) {
        for j in (0..64) {
            *cmap
                .entry(a.chars().nth(i * 64 + j).unwrap())
                .or_insert_with(|| [0; N])
                .get_mut(i)
                .unwrap() |= (1_u64.wrapping_shl(j as u32));
        }
    }

    for i in (0..tlen) {
        // let cmapc = cmap.clone();

        *cmap
            .entry(a.chars().nth(tmax * 64 + i).unwrap())
            .or_insert_with(|| [0; N])
            .get_mut(tmax).unwrap() |= (1_u64.wrapping_shl(i as u32));
            // .unwrap()
            // .get_mut(tmax)
            // .unwrap() |= (1_u64.wrapping_shl(i as u32));
    }

    // log::debug!("tmax={}, tlen={}, a={}, b={}, cmap.keys()={:#?}, asize={}, bsize={}", tmax, tlen, a, b, cmap.keys(), asize, bsize);

    edit_distance_bpv(
        &mut cmap,
        &b.chars().collect::<Vec<char>>(),
        bsize,
        tmax,
        tlen,
    )
}

pub(crate) fn edit_distance<'s>(mut a: &'s str, mut asize: usize, mut b: &'s str, mut bsize: usize) -> usize {
    if asize == 0 {
        return bsize;
    } else if bsize == 0 {
        return asize;
    }

    if asize < bsize {
        (a, b) = (b, a);

        (asize, bsize) = (bsize, asize);
    }

    let mut vsize = ((asize - 1).wrapping_shr(6 as u32)) + 1;
    if vsize > 10 {
        (a, b) = (b, a);
        (asize, bsize) = (bsize, asize);
        vsize = ((asize - 1).wrapping_shr(6 as u32)) + 1;
    }

    match vsize {
        1 => edit_distance_map_::<1>(a, asize, b, bsize),
        2 => edit_distance_map_::<2>(a, asize, b, bsize),
        3 => edit_distance_map_::<3>(a, asize, b, bsize),
        4 => edit_distance_map_::<4>(a, asize, b, bsize),
        5 => edit_distance_map_::<5>(a, asize, b, bsize),
        6 => edit_distance_map_::<6>(a, asize, b, bsize),
        7 => edit_distance_map_::<7>(a, asize, b, bsize),
        8 => edit_distance_map_::<8>(a, asize, b, bsize),
        9 => edit_distance_map_::<9>(a, asize, b, bsize),
        10 => edit_distance_map_::<10>(a, asize, b, bsize),
        _ => edit_distance_dp(a, asize, b, bsize),
    }
}

pub(crate) fn edit_distance_from_str(a: &str, b: &str) -> usize {
    edit_distance(a, a.chars().count(), b, b.chars().count())
}

#[cfg(test)]
mod test {
    use std::time::Instant;

    use super::edit_distance;

    #[test]
    fn t1() {
        println!("{}", 1_i32 << 5 - 1);
    }

    #[test]
    fn t2() {
        println!("{}", !0_u64);
        println!("{}", u64::MAX);
    }

    #[test]
    fn editdistance_test() {
        let str1 = vec![
            "CCTATCAGGGAGCTGTGGGCCAGCCAGGAGGCAGCACATGCCCAATCCCAGGCCCCTCCCGTTGTAAGTTCCCGTTCTACCCGACAGGGACCTGCTGACAAAAGACAGGGCTGGAGAGCCAGCCTGAAGGCCCTGGGACCCTTCTATCCAC",
        "ACTTATGTTTTTAAATGAGGATTATTGATAGTACTCTTGGTTTTTATACCATTCAGATCACTGAATTTATAAAGTACCCATCTAGTACTTCAAAAAGTAAAGTGTTCTGCCAGATCTTAGGTATAGAGGACCCTAACACAGTAAGATCGGA",
        "TAGGGGTATGAGTAGAGCTGAGCTGGGGGAAAAGAGGGAAATTCCCAGGGGTGGAGGAAGAGTCAAGTCCCCCTCTACACCTAGAGGATGAACTTAAGGAAGGAGTGAAGGTCATATGTGTTGTTCCTGAGGAAAAGGCCGCTGTAGAAAA",
        ];

        let str2 = vec![
            "CCTATCAGGGAGCTGTGGGCCAGCCAGGAGGCAGCACATGCCCAATCCCAGGCCCCTCCCGTTGTAAGTTCCCGTTCTACCCGACAGGGACCTGCTGACAAAAGACAGGGCTGGAGAGCCAGCCTGAAGGCCCTGGGACCCTTCTATCCAC",
            "ACTTATGTTTTTAAATGAGGATTATTGATAGTACTCTTGGTTTTTATACCATTCAGATCACTGAATTTATAAAGTACCCATCTAGTACTTGAAAAAGTAAAGTGTTCTGCCAGATCTTAGGTATAGAGGACCCTAACACAGTAAGATCGGA",
            "CCTGGGCCTGGCCCTTGTCTAAAACTGACTCTTTTGAGGGTGATTTTGGATGTTCTTAGTAGAGTCTCTCACCTGTACTTTCCTTGCCTAAGGTGCTGTCTTCTCTTGCAGGTTGCCTACACGTTCCTCACATGCCCTAAGAACCATGGGA",
        ];

        let result = vec![0, 1, 90];

        for i in (0..3) {
            let mut ret = 0;
            let t1 = Instant::now();

            for p in (0..1) {
                ret = edit_distance(
                    str1[i],
                    str1[i].chars().count(),
                    str2[i],
                    str2[i].chars().count(),
                );
            }
            let elapsed = t1.elapsed();
            println!(
                "test 100000 edit_distance, takes {} s\n",
                (elapsed).as_secs_f32()
            );

            if (ret != result[i]) {
                println!(
                    "Fail: (edit_distance), expect {}, but got {}: \n{}\n{}\n",
                    result[i], ret, str1[i], str2[i]
                );
            }
        }
    }

    #[test]
    fn bitwise_reference() {
        println!("{:b}", 3 << 1);
        println!("{:b}", &3 << 2);
    }

    #[test]
    fn unsigned_and_signed_bitwise() {
        println!("{}", (1_i64 >> 3) as u64);
        println!("{}", (1_u64 >> 3) as u64);
        println!("{}", 0_u64 | (1_i64 << 5) as u64);
        println!("{}", 0_u64 | (1_i64 << 5) as u64);
    }

    #[test]
    fn shrrr() {
        println!("{}", (34_i32).wrapping_shr(6));
        println!("{}", (34_i32) >> 6);
    }

}
