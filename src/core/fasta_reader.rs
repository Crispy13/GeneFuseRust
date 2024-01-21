use std::{
    collections::HashMap,
    error::Error,
    fs::File,
    io::{BufRead, BufReader, Read, self},
    mem,
    path::Path,
    fmt,
};

use flate2::read::MultiGzDecoder;

use crate::aux::he::{ErrorExplained, OrExaplain, make_custom_error};

make_custom_error!(EmptyFileError, "Loaded file is empty.");

pub(crate) struct FastaReader {
    m_fasta_file: String,
    m_fasta_file_stream: BufReader<MultiGzDecoder<File>>,
    m_force_upper_case: bool,

    read_buf: Vec<u8>,

    pub(crate) m_current_sequence: String,
    pub(crate) m_current_id: String,
    pub(crate) m_current_description: String,
    pub(crate) m_all_contigs: HashMap<String, String>,
}

impl FastaReader {
    pub(crate) fn new(fasta_file: impl AsRef<Path>, force_upper_case: bool) -> Result<Self, Box<dyn Error>> {
        let fasta_file = fasta_file.as_ref();

        if fasta_file.is_dir() {
            Err(format!(
                "There is a problem with the provided fasta file: \
            '{}' is a directory NOT a file...\n",
                fasta_file.to_str().unwrap()
            ))?
        }

        let mut fasta_buf_reader = BufReader::new(MultiGzDecoder::new(File::open(&fasta_file)?));

        // TODO:verify that the file can be read

        // seek to first contig
        let mut read_buf = Vec::new();
        if let Ok(true) = fasta_buf_reader
            .read_until(b'>', &mut read_buf)
            .and_then(|rl| Ok(rl > 0))
        {
            // do nothing
        } else {
            Err(EmptyFileError::new())?
        }

        read_buf.clear();
        
        Ok(Self {
            m_fasta_file: fasta_file.to_str().unwrap().to_string(),
            m_fasta_file_stream: fasta_buf_reader,
            m_force_upper_case: force_upper_case,
            m_current_sequence: String::new(),
            m_current_id: String::new(),
            m_current_description: String::new(),
            m_all_contigs: HashMap::new(),

            read_buf: read_buf,
        })
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

    fn contigs(&self) -> &HashMap<String, String> {
        &self.m_all_contigs
    }

    fn read_next(&mut self) -> bool {
        const header_delimiter: [u8; 2] = [b'\n', b' '];

        let read_buf = &mut self.read_buf;
        let mut found_header = false;
        let m_force_upper_case = &mut self.m_force_upper_case;

        let mut ss_header = String::new();
        let mut ss_seq = String::new();

        let mut has_read_somethings = false;
        // read_buf.clear();
        // read_buf.push(b'>'); // we stopped at first '>' in `new()`. so add it to the front first.
        if let Ok(true) = self
            .m_fasta_file_stream
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

        self.m_current_id = ss_header;
        self.m_current_sequence = ss_seq;

        has_read_somethings
    }

    pub(crate) fn read_all(&mut self) {
        while self.read_next() {
            self.m_all_contigs.insert(
                mem::take(&mut self.m_current_id),
                mem::take(&mut self.m_current_sequence),
            );
        }
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
        reader.m_fasta_file_stream.read_to_string(&mut s).unwrap();

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
