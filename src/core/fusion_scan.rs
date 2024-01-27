use std::error::Error;

use super::{fusion::Fusion, pescanner::PairEndScanner};

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

    pub(crate) fn scan(self) -> Result<bool, Box<dyn Error>> {
        // let fusions = Fusion::parse_csv(&self.m_fusion_file)?;

        if self.m_read2_file != "" {
            let mut pescanner = PairEndScanner::new(
                self.m_fusion_file,
                self.m_ref_file,
                self.m_read1_file,
                self.m_read2_file,
                self.m_html_file,
                self.m_json_file,
                self.m_thread_num as i32,
            );
            
            Ok(pescanner.scan()?)
        } else {
            unimplemented!()
            // SingleEndScanner sescanner( mFusionFile, mRefFile, mRead1File, mHtmlFile, mJsonFile, mThreadNum);
            // return sescanner.scan();
        }

    }
}
