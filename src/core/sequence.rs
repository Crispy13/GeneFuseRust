#[derive(Debug, PartialEq, PartialOrd, Clone)]
pub(crate) struct Sequence {
    pub(crate) m_str: String,
}

impl Sequence {
    pub(crate) fn new(m_str: String) -> Self {
        Self {
            m_str,
        }
    }

    pub(crate) fn reverse_complement(&self) -> Self {
        Self { m_str: reverse_complement(&self.m_str) }
    }

    pub(crate) fn len(&self) -> usize {
        self.m_str.len()
    }
}

pub(crate) fn reverse_complement(seq:&str) -> String {
    let mut s = seq.chars().collect::<Vec<char>>();
    let mut s_iter = s.iter_mut();

    // let (mut f, mut b);

    loop {
        match (s_iter.next(), s_iter.next_back()) {
            (Some(f), Some(b)) => {
                (*f, *b) = (get_complement_base(b), get_complement_base(f));
            },
            (Some(f), None) => {
                *f = get_complement_base(f);
                break;
            },
            (None, Some(b)) => { // this is unreachable.
                unreachable!();
                *b = get_complement_base(b);
                break;
            },
            (None, None) => {
                break;
            }
        }
    }

    s.into_iter().collect::<String>()

}

fn get_complement_base(base: &char) -> char {
    match base {
        'A' | 'a' => 'T',
        'T' | 't' => 'A',
        'C' | 'c' => 'G',
        'G' | 'g' => 'C',
        _=> 'N',
    }
}

#[cfg(test)]
mod test {
    use super::{get_complement_base, Sequence};

    #[test]
    fn re_test() {
        assert_eq!("AACCCGCAT", Sequence::new("ATGCGGGTT".into()).reverse_complement().m_str.as_str());
        assert_eq!("CTANTTCG", Sequence::new("CGAANTAG".into()).reverse_complement().m_str.as_str());
    }
}