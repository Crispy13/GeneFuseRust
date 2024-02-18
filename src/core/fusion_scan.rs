use std::{
    error,
    fs::File,
    io::{BufRead, BufReader},
    marker::PhantomData,
    mem,
    path::{Path, PathBuf},
    process::exit,
    sync::Arc,
};

use rayon::{prelude::*, ThreadPoolBuilder};

use crate::{
    aux::limited_bufreader::LimitedBufReader,
    core::{
        fasta_reader::FastaReader, fastq_reader::{FastqReader, FastqReaderPair}, read::{SequenceReadCow, SequenceReadPairCow}, sescanner::{self, SingleEndScanner}
    },
};

use super::{fusion::Fusion, pescanner::PairEndScanner};

pub(crate) type Error = Box<dyn error::Error + Send + Sync>;
pub(crate) struct FusionScan {
    m_fusion_file: String,
    m_read1_file: String,
    m_read2_file: String,
    m_html_file: String,
    m_json_file: String,
    m_ref_file: String,
    m_thread_num: usize,
}

impl FusionScan {
    pub(crate) fn new(
        m_fusion_file: String,
        m_ref_file: String,
        m_read1_file: String,
        m_read2_file: String,
        m_html_file: String,
        m_json_file: String,

        m_thread_num: usize,
    ) -> FusionScan {
        Self {
            m_fusion_file,
            m_read1_file,
            m_read2_file,
            m_html_file,
            m_json_file,
            m_ref_file,
            m_thread_num,
        }
    }

    pub(crate) fn scan_per_fusion_csv(&self) -> Result<bool, Error> {
        // read reference first.
        let ref_file = self.m_ref_file.as_str();
        let mut m_reference = FastaReader::new(ref_file, false)?;

        log::debug!("Reading reference, {}", &ref_file);
        m_reference.read_all();


        // read input seq fastqs
        log::info!("Reading input seqeunces...");
        let (srp_vec, sr_vec) = if !self.m_read2_file.is_empty() {
            let mut fqr = FastqReaderPair::from_paths(&self.m_read1_file, &self.m_read2_file)?;
            
            let mut srp_vec = vec![];
            while let Some(srp) = fqr.read() {
                srp_vec.push(srp);
            }

            (Some(srp_vec), None)
        } else {
            let mut fqr = FastqReader::new(&self.m_read1_file, true)?;
            
            let mut sr_vec = vec![];
            while let Some(srp) = fqr.read() {
                sr_vec.push(srp);
            }

            (None, Some(sr_vec))
        };

        // if let Some(s) = srp_vec {
        //     let a = s.iter().cloned().map(SequenceReadPairCow::Owned).collect::<Vec<_>>();
        //     let b = s.iter().map(SequenceReadPairCow::Borrowed).collect::<Vec<_>>();

        //     let c = b.iter().map(|e| {
        //         (SequenceReadCow::Owned(e.m_left.clone()), SequenceReadCow::Owned(e.m_right.clone()))
        //     }).collect::<Vec<_>>();

        //     exit(0);
        // }

        let mut scanner_m_ref = ScannerFastaReader::new(m_reference);

        let fusion_csv_paths = self.get_fusion_csv_vec_from_input()?;
        let (html_file_paths, json_file_paths) =
            self.get_report_names_from_fusion_csvs(fusion_csv_paths.as_slice());

        let (outer_thread_num, inner_thread_num) = if fusion_csv_paths.len() >= self.m_thread_num {
            (self.m_thread_num, 1_usize)
        } else {
            (
                fusion_csv_paths.len(),
                self.m_thread_num / fusion_csv_paths.len(),
            )
        };

        log::info!("given csv count={}, parallel job count={}, inner_thread_num={}", fusion_csv_paths.len(), outer_thread_num, inner_thread_num);

        let tp = ThreadPoolBuilder::new()
            .num_threads(outer_thread_num)
            .thread_name(|i| format!("MainThreadPool-{i}"))
            .build()
            .unwrap();

        tp.install(|| {
            fusion_csv_paths
                .into_iter()
                .zip(html_file_paths.into_iter().zip(json_file_paths))
                .par_bridge()
                .map(|(fusion_csv, (html_file, json_file))| {
                    let res = if self.m_read2_file != "" {
                        let pescanner = PairEndScanner::new(
                            fusion_csv,
                            self.m_ref_file.clone(),
                            self.m_read1_file.clone(),
                            self.m_read2_file.clone(),
                            html_file,
                            json_file,
                            inner_thread_num as i32,
                            srp_vec.as_ref().map(|v| v.as_slice()),
                        );

                        scanner_m_ref.scan_per_fusion_csv(pescanner)
                    } else {
                        let sescanner = SingleEndScanner::new(
                            fusion_csv,
                            self.m_ref_file.clone(),
                            self.m_read1_file.clone(),
                            html_file,
                            json_file,
                            inner_thread_num as i32,
                            sr_vec.as_ref().map(|v| v.as_slice()),
                        );

                        scanner_m_ref.scan_per_fusion_csv(sescanner)
                    };

                    res
                })
                .collect::<Result<Vec<bool>, Error>>()
                .and_then(|vb| Ok(vb.into_iter().all(|e| e)))
        })
    }

    fn get_report_names_from_fusion_csvs(
        &self,
        fusion_csv_paths: &[String],
    ) -> (Vec<String>, Vec<String>) {
        let hf = self.m_html_file.as_str();

        let (hf_parent, hf_stem, hf_ext) = if !hf.is_empty() {
            let hf = Path::new(hf);
            let (hf_parent, hf_stem, hf_ext) = (
                hf.parent().unwrap().to_str().unwrap(),
                hf.file_stem().unwrap().to_str().unwrap(),
                hf.extension().unwrap().to_str().unwrap(),
            );
            (hf_parent, hf_stem, hf_ext)
        } else {
            ("", "", "")
        };

        let jf = self.m_json_file.as_str();

        let (jf_parent, jf_stem, jf_ext) = if !jf.is_empty() {
            let jf = Path::new(jf);
            let (jf_parent, jf_stem, jf_ext) = (
                jf.parent().unwrap().to_str().unwrap(),
                jf.file_stem().unwrap().to_str().unwrap(),
                jf.extension().unwrap().to_str().unwrap(),
            );
            (jf_parent, jf_stem, jf_ext)
        } else {
            ("", "", "")
        };

        let mut html_file_vec = Vec::new();
        let mut json_file_vec = Vec::new();

        fusion_csv_paths.iter().for_each(|fc| {
            let fc_stem = Path::new(fc).file_stem().unwrap().to_str().unwrap();
            if !hf.is_empty() {
                html_file_vec.push(
                    [hf_parent, &format!("{}_{}.{}", hf_stem, fc_stem, hf_ext)]
                        .iter()
                        .collect::<PathBuf>()
                        .into_os_string()
                        .into_string()
                        .unwrap(),
                );
            }

            if !jf.is_empty() {
                json_file_vec.push(
                    [jf_parent, &format!("{}_{}.{}", jf_stem, fc_stem, jf_ext)]
                        .iter()
                        .collect::<PathBuf>()
                        .into_os_string()
                        .into_string()
                        .unwrap(),
                );
            }
        });

        (html_file_vec, json_file_vec)
    }

    fn get_fusion_csv_vec_from_input(&self) -> Result<Vec<String>, Error> {
        let mut f = LimitedBufReader::new(
            BufReader::new(File::open(self.m_fusion_file.as_str())?),
            1000,
        );

        let mut fusion_csvs = Vec::new();

        let mut buf_s = String::new();
        while let Ok(true) = {
            buf_s.clear();
            f.read_line(&mut buf_s).and_then(|rl| Ok(rl > 0))
        } {
            let s = buf_s.trim();
            if s.is_empty() {
                continue;
            }

            if !Path::new(s).is_file() {
                eprintln!("Fusion csv file '{s}' was not found.");
                exit(-1);
            } else {
                fusion_csvs.push(s.to_string());
            }
        }

        Ok(fusion_csvs)
    }

    fn scan_single_csv(self) -> Result<bool, Error> {
        if self.m_read2_file != "" {
            let mut pescanner = PairEndScanner::new(
                self.m_fusion_file,
                self.m_ref_file,
                self.m_read1_file,
                self.m_read2_file,
                self.m_html_file,
                self.m_json_file,
                self.m_thread_num as i32,
                None,
            );

            Ok(pescanner.scan()?)
        } else {
            let mut sescanner = SingleEndScanner::new(
                self.m_fusion_file,
                self.m_ref_file,
                self.m_read1_file,
                self.m_html_file,
                self.m_json_file,
                self.m_thread_num as i32,
                None,
            );

            Ok(sescanner.scan()?)
        }
    }

    pub(crate) fn scan(self) -> Result<bool, Error> {
        // run proper function by the type of input fusion file.
        match Path::new(self.m_fusion_file.as_str())
            .extension()
            .unwrap()
            .to_str()
            .unwrap()
        {
            "csv" => self.scan_single_csv(),
            _ => self.scan_per_fusion_csv(),
        }
    }
}

struct ScannerFastaReader {
    fasta_reader: Option<Arc<FastaReader>>, // use Option to use mem::take.
}

impl ScannerFastaReader {
    fn new(fasta_reader: FastaReader) -> Self {
        Self {
            fasta_reader: Some(Arc::new(fasta_reader)),
        }
    }
}

impl ScannerFastaReader {
    fn scan_per_fusion_csv<S: Scanner>(&self, mut scanner: S) -> Result<bool, Error> {
        // let r = scanner.scan_per_fusion_csv(mem::take(&mut self.fasta_reader).unwrap());

        let r = scanner.scan_per_fusion_csv(Arc::clone(self.fasta_reader.as_ref().unwrap()));

        // self.fasta_reader = Some(scanner.drop_and_get_back_fasta_reader());

        r
    }
}

pub(crate) trait Scanner {
    fn scan(&mut self) -> Result<bool, Error>;

    fn scan_per_fusion_csv(&mut self, ref_fasta: Arc<FastaReader>) -> Result<bool, Error>;

    // fn drop_and_get_back_fasta_reader(self) -> FastaReader;
}

impl<'s> Scanner for PairEndScanner<'s> {
    fn scan(&mut self) -> Result<bool, Error> {
        self.scan()
    }

    fn scan_per_fusion_csv(&mut self, ref_fasta: Arc<FastaReader>) -> Result<bool, Error> {
        self.scan_per_fusion_csv(ref_fasta)
    }

    // fn drop_and_get_back_fasta_reader(self) -> FastaReader {
    //     self.drop_and_get_back_fasta_reader()
    // }
}

impl<'s> Scanner for SingleEndScanner<'s> {
    fn scan(&mut self) -> Result<bool, Error> {
        self.scan()
    }

    fn scan_per_fusion_csv(&mut self, ref_fasta: Arc<FastaReader>) -> Result<bool, Error> {
        self.scan_per_fusion_csv(ref_fasta)
    }

    // fn drop_and_get_back_fasta_reader(self) -> FastaReader {
    //     self.drop_and_get_back_fasta_reader()
    // }
}
