pub(crate) mod logging;

use std::{path::Path, process::exit};

#[inline]
pub(crate) fn check_file_valid(s: impl AsRef<Path>) -> () {
    let s = s.as_ref();

    if !s.is_file() {
        println!(
            "ERROR: file '{}' doesn't exist, quit now",
            s.to_str().unwrap()
        );
        exit(-1);
    }

    if s.is_dir() {
        println!(
            "ERROR: '{}' is a folder, not a file, quit now",
            s.to_str().unwrap()
        );
        exit(-1);
    }
}

pub(crate) trait StringCPP {
    fn subchars(&self, pos: usize, n: usize) -> &str;
}

impl StringCPP for String {
    fn subchars(&self, pos: usize, n: usize) -> &str {
        // self.chars().skip(pos).take(n).collect::<String>()
        self.get(pos..(pos+n)).unwrap()
    }
}

impl StringCPP for &str {
    fn subchars(&self, pos: usize, n: usize) -> &str {
        self.get(pos..(pos+n)).unwrap()
    }
}
#[inline]
pub(crate) fn dis_connected_count(s:&str) -> i32 {
    let mut diff = 0;
    for i in (0..(s.len()-1)) {
        if s.get(i..(i+1)).unwrap() != s.get((i+1)..(i+2)).unwrap() {
            diff +=1;
        }
    }
    return diff
}

#[inline]
pub(crate) fn int2str(num: i32) -> String {
    num.to_string()
}


