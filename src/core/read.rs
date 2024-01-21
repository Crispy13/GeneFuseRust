use super::sequence::Sequence;

#[derive(Debug, PartialEq, PartialOrd)]
pub(crate) struct SequenceRead {
    pub(crate) m_name: String,
    pub(crate) m_seq: Sequence,
    pub(crate) m_strand: String,
    pub(crate) m_quality: String,
    pub(crate) m_has_quality: bool,
}

impl SequenceRead {
    pub(crate) fn new(
        m_name: String,
        m_seq: String,
        m_strand: String,
        m_quality: String,
        m_has_quality: bool,
    ) -> Self {
        Self {
            m_name,
            m_seq: Sequence::new(m_seq),
            m_strand,
            m_quality,
            m_has_quality,
        }
    }

    pub(crate) fn len(&self) {
        self.m_seq.len()
    }
}

pub(crate) struct SequenceReadPair {
    m_left: SequenceRead,
    m_right: SequenceRead,
}

impl SequenceReadPair {
    pub(crate) fn new(left: SequenceRead, right: SequenceRead) -> SequenceReadPair {
        Self {
            m_left: left,
            m_right: right,
        }
    }
}
