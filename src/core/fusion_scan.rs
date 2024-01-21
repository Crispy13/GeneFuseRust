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

    fn scan(&self) {
        // let fusions = Fusion::parse_csv(&self.m_fusion_file);
        // if(mRead2File != ""){
        //     PairEndScanner pescanner( mFusionFile, mRefFile, mRead1File, mRead2File, mHtmlFile, mJsonFile, mThreadNum);
        //     return pescanner.scan();
        // }
        // else{
        //     SingleEndScanner sescanner( mFusionFile, mRefFile, mRead1File, mHtmlFile, mJsonFile, mThreadNum);
        //     return sescanner.scan();
        // }
    }
}
