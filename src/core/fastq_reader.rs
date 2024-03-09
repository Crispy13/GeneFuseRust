use super::{fusion_scan::Error, read::SequenceReadPairArc};

use std::{
    borrow::Cow,
    error,
    fs::File,
    io::{BufRead, BufReader, Read},
    path::Path,
    process::exit,
};

use flate2::read::MultiGzDecoder;
use log4rs::append::file;

use crate::aux::limited_bufreader::LimitedBufReader;

use super::read::{SequenceRead, SequenceReadPair};

pub(crate) struct FastqReader {
    m_filename: String,
    m_zip_file: Option<LimitedBufReader<BufReader<MultiGzDecoder<File>>>>,
    m_file: Option<LimitedBufReader<BufReader<File>>>,
    m_zipped: bool,
    m_has_quality: bool,
}

pub(crate) const max_take: u64 = 1000;

impl FastqReader
// where
//     R: BufRead,
{
    pub(crate) fn new(
        file_name: impl AsRef<Path>,
        has_quality: bool,
    ) -> Result<Self, Error> {
        let file_name = file_name.as_ref();

        let (mut gzr, mut nr) = (None, None);
        let m_zipped;
        if Self::is_zip_fastq(file_name) {
            m_zipped = true;
            gzr = Some(LimitedBufReader::new(
                BufReader::new(MultiGzDecoder::new(File::open(file_name)?)),
                max_take,
            ));
        } else if Self::is_fastq(file_name) {
            m_zipped = false;
            nr = Some(LimitedBufReader::new(
                BufReader::new(File::open(file_name)?),
                max_take,
            ));
        } else {
            eprintln!("ERROR: the input file should be fastq (.fq, .fastq) or gzipped fastq (.fq.gz, .fastq.gz) {}", file_name.to_str().unwrap());
            exit(-1);
        };

        Ok(Self {
            m_filename: file_name.to_str().unwrap().to_string(),
            m_zip_file: gzr,
            m_file: nr,
            m_zipped,
            m_has_quality: has_quality,
        })
    }

    fn buf_reader(&mut self) -> &mut dyn BufRead {
        if let Some(gz_buf_reader) = self.m_zip_file.as_mut() {
            gz_buf_reader
        } else {
            self.m_file.as_mut().unwrap()
        }
    }

    pub(crate) fn read(&mut self) -> Option<SequenceRead> {
        let m_has_quality = self.m_has_quality;
        let buf_reader = self.buf_reader();

        let mut s = String::new();
        let name = if let Ok(true) = buf_reader.read_line(&mut s).and_then(|rl| Ok(rl > 0)) {
            let mut c = s.clone();
            c.pop().and_then(|ch| { // remove newline character if it exists
                if ch != '\n' {
                    c.push(ch);
                }
                Some(())
            });
            s.clear();
            c
        } else {
            return None;
        };

        let sequence = if let Ok(true) = buf_reader.read_line(&mut s).and_then(|rl| Ok(rl > 0)) {
            let mut c = s.clone();
            c.pop().and_then(|ch| {
                if ch != '\n' {
                    c.push(ch);
                }
                Some(())
            });
            s.clear();
            c
        } else {
            return None;
        };

        let strand = if let Ok(true) = buf_reader.read_line(&mut s).and_then(|rl| Ok(rl > 0)) {
            let mut c = s.clone();
            c.pop().and_then(|ch| {
                if ch != '\n' {
                    c.push(ch);
                }
                Some(())
            });
            s.clear();
            c
        } else {
            return None;
        };

        let quality = if m_has_quality {
            if let Ok(true) = buf_reader.read_line(&mut s).and_then(|rl| Ok(rl > 0)) {
                let mut c = s.clone();
                c.pop().and_then(|ch| {
                    if ch != '\n' {
                        c.push(ch);
                    }
                    Some(())
                });
                s.clear();
                c
            } else {
                return None;
            }
        } else {
            String::new()
        };

        Some(SequenceRead::new(
            name,
            sequence,
            strand,
            quality,
            m_has_quality,
        ))
    }

    fn is_zip_fastq(file_name: &Path) -> bool {
        let file_name = file_name.to_str().unwrap();

        if file_name.ends_with(".fastq.gz") {
            return true;
        } else if file_name.ends_with(".fq.gz") {
            return true;
        } else if file_name.ends_with(".fasta.gz") {
            return true;
        } else if file_name.ends_with(".fa.gz") {
            return true;
        } else {
            return false;
        }
    }

    fn is_fastq(file_name: &Path) -> bool {
        let file_name = file_name.to_str().unwrap();

        if file_name.ends_with(".fastq") {
            return true;
        } else if file_name.ends_with(".fq") {
            return true;
        } else if file_name.ends_with(".fasta") {
            return true;
        } else if file_name.ends_with(".fa") {
            return true;
        } else {
            return false;
        }
    }

    fn is_zipped(&self) -> bool {
        self.m_zipped
    }
}

pub(crate) struct FastqReaderPair {
    m_left: FastqReader,
    m_right: FastqReader,
}

impl FastqReaderPair {
    pub(crate) fn new(left: FastqReader, right: FastqReader) -> FastqReaderPair {
        Self {
            m_left: left,
            m_right: right,
        }
    }

    pub(crate) fn from_paths(
        left_name: impl AsRef<Path>,
        right_name: impl AsRef<Path>,
    ) -> Result<FastqReaderPair, Error> {
        Ok(Self {
            m_left: FastqReader::new(left_name, true)?,
            m_right: FastqReader::new(right_name, true)?,
        })
    }

    pub(crate) fn read(&mut self) -> Option<SequenceReadPair> {
        let l = self.m_left.read();
        let r = self.m_right.read();

        if l.is_none() || r.is_none() {
            None
        } else {
            Some(SequenceReadPair::new(l.unwrap(), r.unwrap()))
        }
    }

    pub(crate) fn read_and_into_arc(&mut self) -> Option<SequenceReadPairArc> {
        self.read().and_then(|sp| Some(SequenceReadPairArc::from_seq_read_pair(sp)))
    }
}

#[cfg(test)]
mod test {
    use std::path::PathBuf;

    use super::*;

    struct A {
        gzr: BufReader<MultiGzDecoder<File>>,
        nr: BufReader<File>,
    }

    fn _var_init() -> Result<(), Error> {
        let file_name = PathBuf::from("testdata/R1.fq");
        let (gzr, nr);
        let m_zipped;
        let buf_reader: &dyn BufRead = if FastqReader::is_zip_fastq(&file_name) {
            m_zipped = true;
            gzr = BufReader::new(MultiGzDecoder::new(File::open(file_name)?));
            &gzr
        } else {
            m_zipped = false;
            nr = BufReader::new(File::open(file_name)?);
            &nr
        };

        // println!("{:?} {:?}", gzr, nr);

        // A {
        //     gzr,
        //     nr,
        // };

        Ok(())
    }

    #[test]
    fn var_init() {
        _var_init().unwrap();
    }

    #[test]
    fn new_fq_reader() {
        FastqReader::new("testdata/R1.fq", true).unwrap();
    }

    #[test]
    fn fq_test() {
        println!("{}", _test());
    }

    fn _test() -> bool {
        let mut reader1 = FastqReader::new("testdata/R1.fq", true).unwrap();
        let mut reader2 = FastqReader::new("testdata/R1.fq.gz", true).unwrap();

        let (mut r1, mut r2);
        loop {
            r1 = reader1.read();
            r2 = reader2.read();

            if r1.is_none() || r2.is_none() {
                break;
            }

            println!("{:?}", r1.as_ref().unwrap());
            println!("{:?}", r2.as_ref().unwrap());

            if &r1.unwrap().m_seq.m_str != &r2.unwrap().m_seq.m_str {
                return false;
            }
        }

        return true;
    }
}
