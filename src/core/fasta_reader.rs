use std::panic::Location;
use std::{
    collections::{BTreeMap, HashMap},
    error,
    fmt,
    fs::File,
    io::{self, BufRead, BufReader, Read},
    mem,
    path::Path,
};

use flate2::read::MultiGzDecoder;

use crate::aux::he::{
    make_custom_error, make_custom_error3, make_custom_error4, ErrorExplained, OrExaplain,
};
use crate::aux::limited_bufreader::LimitedBufReader;
use crate::aux::pbar::prepare_pbar;

use super::fusion_scan::Error;

make_custom_error4!(EmptyFileError, "Loaded file is empty.");

pub(crate) struct FastaReader {
    pub(crate) m_fasta_file: String,
    // m_fasta_file_stream: BufReader<MultiGzDecoder<File>>,
    m_fasta_gz_reader: Option<BufReader<MultiGzDecoder<File>>>,
    m_fasta_txt_reader: Option<BufReader<File>>,
    m_force_upper_case: bool,

    read_buf: Vec<u8>,

    pub(crate) m_current_sequence: String,
    pub(crate) m_current_id: String,
    pub(crate) m_current_description: String,
    pub(crate) m_all_contigs: BTreeMap<String, String>,
}

impl FastaReader {
    pub(crate) fn new(
        fasta_file: impl AsRef<Path>,
        force_upper_case: bool,
    ) -> Result<Self, Error> {
        let fasta_file = fasta_file.as_ref();

        if fasta_file.is_dir() {
            Err(format!(
                "There is a problem with the provided fasta file: \
            '{}' is a directory NOT a file...\n",
                fasta_file.to_str().unwrap()
            ))?
        }

        // const max_take:u64 = 10000000;
        let (mut gz_reader, mut txt_reader) = (None, None);
        // load fasta
        let fasta_buf_reader: &mut dyn BufRead = if fasta_file.extension().unwrap() == "gz" {
            // gzipped fasta
            gz_reader = Some(BufReader::new(MultiGzDecoder::new(File::open(
                &fasta_file,
            )?)));

            gz_reader.as_mut().unwrap()
        } else {
            // text fasta
            txt_reader = Some(BufReader::new(File::open(&fasta_file)?));
            txt_reader.as_mut().unwrap()
        };

        // let mut fasta_buf_reader = BufReader::new(MultiGzDecoder::new(File::open(&fasta_file)?));

        // TODO:verify that the file can be read

        // seek to first contig
        let mut read_buf = Vec::new();
        if let Ok(true) = fasta_buf_reader
            .read_until(b'>', &mut read_buf)
            .and_then(|rl| Ok(rl > 0))
        {
            // do nothing
        } else {
            Err(EmptyFileError::new(&fasta_file))?
        }

        read_buf.clear();

        Ok(Self {
            m_fasta_file: fasta_file.to_str().unwrap().to_string(),
            m_fasta_gz_reader: gz_reader,
            m_fasta_txt_reader: txt_reader,
            m_force_upper_case: force_upper_case,
            m_current_sequence: String::new(),
            m_current_id: String::new(),
            m_current_description: String::new(),
            m_all_contigs: BTreeMap::new(),

            read_buf: read_buf,
        })
    }

    fn m_fasta_file_stream(&mut self) -> &mut dyn BufRead {
        if let Some(gzr) = self.m_fasta_gz_reader.as_mut() {
            gzr
        } else {
            self.m_fasta_txt_reader.as_mut().unwrap()
        }
    }

    fn current_id(&self) -> &str {
        self.m_current_id.as_str()
    }

    fn current_description(&self) -> &str {
        self.m_current_description.as_str()
    }

    fn current_sequence(&self) -> &str {
        self.m_current_sequence.as_str()
    }

    fn contigs(&self) -> &BTreeMap<String, String> {
        &self.m_all_contigs
    }

    fn read_next(&mut self) -> bool {
        const header_delimiter: [u8; 2] = [b'\n', b' '];

        let read_buf = &mut self.read_buf;
        let mut found_header = false;
        let m_force_upper_case = &mut self.m_force_upper_case;

        let fasta_buf_reader: &mut dyn BufRead = if self.m_fasta_gz_reader.is_some() {
            self.m_fasta_gz_reader.as_mut().unwrap()
        } else {
            self.m_fasta_txt_reader.as_mut().unwrap()
        };

        let mut ss_header = String::new();
        let mut ss_seq = String::new();

        let mut has_read_somethings = false;
        // read_buf.clear();
        // read_buf.push(b'>'); // we stopped at first '>' in `new()`. so add it to the front first.
        if let Ok(true) = fasta_buf_reader
            .read_until(b'>', read_buf)
            .and_then(|rl| Ok(rl > 0))
        {
            read_buf.pop().and_then(|b| {
                if b != b'>' {
                    read_buf.push(b);
                }

                Some(())
            });

            has_read_somethings = true;
            let mut read_buf_drained = read_buf.drain(..);

            while let Some((true, b)) = read_buf_drained
                .next()
                .and_then(|b| Some((!header_delimiter.contains(&b), b)))
            {
                ss_header.push(char::from(b));
            }

            if *m_force_upper_case {
                ss_seq.extend(read_buf_drained.into_iter().filter_map(|b| {
                    if b != b'\n' {
                        filter_map_valid_seq_to_upper(b)
                    } else {
                        None
                    }
                }))
            } else {
                ss_seq.extend(read_buf_drained.into_iter().filter_map(|b| {
                    if b != b'\n' {
                        filter_map_valid_seq(b)
                    } else {
                        None
                    }
                }))
            }
        }

        // log::debug!("read header = {}, current_seq = {} ..", ss_header, ss_seq.chars().take(20).collect::<String>());

        self.m_current_id = ss_header;
        self.m_current_sequence = ss_seq;

        has_read_somethings
    }

    pub(crate) fn read_all(&mut self) {
        let pbar = prepare_pbar(0);
        pbar.set_message("Reading references...");
        while self.read_next() {
            pbar.inc(1);
            self.m_all_contigs.insert(
                mem::take(&mut self.m_current_id),
                mem::take(&mut self.m_current_sequence),
            );
        }

        pbar.finish_and_clear();
    }
}

#[inline]
fn filter_map_valid_seq_to_upper(mut b: u8) -> Option<char> {
    if b.is_ascii_alphabetic() || b == b'-' || b == b'*' {
        if b.is_ascii_lowercase() {
            b.make_ascii_uppercase()
        }
        Some(char::from(b))
    } else {
        None
    }
}

#[inline]
fn filter_map_valid_seq(mut b: u8) -> Option<char> {
    if b.is_ascii_alphabetic() || b == b'-' || b == b'*' {
        Some(char::from(b))
    } else {
        None
    }
}

#[cfg(test)]
mod test {
    use std::{error::Error, io::Read};

    use super::FastaReader;

    // #[test]
    fn _fasta_reader() -> Result<(), Box<dyn Error>> {
        let mut reader = FastaReader::new("testdata/tinyref.fa.gz", true).unwrap();

        reader.read_all();

        let contig1 = "GATCACAGGTCTATCACCCTATTAATTGGTATTTTCGTCTGGGGGGTGTGGAGCCGGAGCACCCTATGTCGCAGT";
        let contig2 = "GTCTGCACAGCCGCTTTCCACACAGAACCCCCCCCTCCCCCCGCTTCTGGCAAACCCCAAAAACAAAGAACCCTA";

        if !reader.m_all_contigs.contains_key("contig1")
            || !reader.m_all_contigs.contains_key("contig2")
        {
            println!("{:#?}", reader.m_all_contigs);
            Err("1")?;
        }

        if reader.m_all_contigs.get("contig1").unwrap() != contig1
            || reader.m_all_contigs.get("contig2").unwrap() != contig2
        {
            println!("{:#?}", reader.m_all_contigs);
            Err("2")?;
        }

        Ok(())
    }

    #[test]
    fn read_next() -> () {
        let mut reader = FastaReader::new("testdata/tinyref.fa.gz", true).unwrap();

        println!("{}", reader.read_next());

        println!("{:#?}", reader.m_all_contigs);
    }

    #[test]
    fn read_to_end() {
        let mut reader = FastaReader::new("testdata/tinyref.fa.gz", true).unwrap();

        let mut s = String::new();
        reader.m_fasta_file_stream().read_to_string(&mut s).unwrap();

        println!("s={s:#}");
    }

    #[test]
    fn fasta_reader() {
        _fasta_reader().unwrap();
    }

    // FastaReader reader("testdata/tinyref.fa");
    // reader.readAll();

    // string contig1 = "GATCACAGGTCTATCACCCTATTAATTGGTATTTTCGTCTGGGGGGTGTGGAGCCGGAGCACCCTATGTCGCAGT";
    // string contig2 = "GTCTGCACAGCCGCTTTCCACACAGAACCCCCCCCTCCCCCCGCTTCTGGCAAACCCCAAAAACAAAGAACCCTA";

    // if(reader.mAllContigs.count("contig1") == 0 || reader.mAllContigs.count("contig2") == 0 )
    //     return false;

    // if(reader.mAllContigs["contig1"] != contig1 || reader.mAllContigs["contig2"] != contig2 )
    //     return false;

    // return true;
}
