pub(crate) mod logging;

use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
    process::exit,
};

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
        self.get(pos..(pos + n)).unwrap()
    }
}

impl StringCPP for &str {
    fn subchars(&self, pos: usize, n: usize) -> &str {
        self.get(pos..(pos + n)).unwrap()
    }
}
#[inline]
pub(crate) fn dis_connected_count(s: &str) -> i32 {
    let mut diff = 0;
    for i in (0..(s.len() - 1)) {
        if s.get(i..(i + 1)).unwrap() != s.get((i + 1)..(i + 2)).unwrap() {
            diff += 1;
        }
    }
    return diff;
}

#[inline]
pub(crate) fn int2str(num: i32) -> String {
    num.to_string()
}

#[inline]
pub(crate) fn open_csv() -> File {
    let file = File::create("object.tsv").unwrap();

    file
}

pub(crate) struct ObjectTsvWriter {
    buf_writer: BufWriter<File>,
}

impl Default for ObjectTsvWriter {
    fn default() -> Self {
        Self {
            buf_writer: BufWriter::new(File::create("object.tsv").unwrap()),
        }
    }
}

impl Write for ObjectTsvWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.buf_writer.write(buf)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.buf_writer.flush()
    }
}

macro_rules! write_tsv_row {
    ($dst:expr, $($x:expr),*) => {
        {
            let mut s = String::new();
            $(
                write!(&mut s, "{}\t", $x).unwrap();
            )*
            
            s.pop().unwrap(); // remove the last ','

            writeln!($dst, "{}", s).unwrap();
    }
    };
}
pub(crate) use write_tsv_row;

#[cfg(test)]
mod test {
    use std::io::Write;
    use std::fmt::Write as fmt_write;

    use super::ObjectTsvWriter;

    #[test]
    fn wcr() {
        let mut a = ObjectTsvWriter::default();

        let b= 2;
        let c= "SDF";
        let f = 2.2;
        write_tsv_row!(&mut a, b,c,f);

    }
}
