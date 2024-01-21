pub(crate) mod logging;

use std::{path::Path, process::exit};

#[inline]
pub(crate) fn check_file_valid(s: impl AsRef<Path>) -> () {
    let s = s.as_ref();

    if !s.is_file() {
        println!("ERROR: file '{}' doesn't exist, quit now", s.to_str().unwrap());
        exit(-1);
    }

    if s.is_dir() {
        println!("ERROR: '{}' is a folder, not a file, quit now", s.to_str().unwrap());
        exit(-1);
    }

}