use core::num;
use std::{error, fmt::Write};

use crate::aux::he::OrExaplain;

use super::fusion_scan::Error;
// use anyhow::Result;

#[derive(Debug, Clone)]
pub(crate) struct Exon {
    id: i32,
    start: i32,
    end: i32,
}

#[derive(Debug, Clone)]
pub(crate) struct Gene {
    pub(crate) m_name: String,
    pub(crate) m_chr: String,
    pub(crate) m_start: i32,
    pub(crate) m_end: i32,
    pub(crate) m_exons: Vec<Exon>,
    pub(crate) m_reversed: bool,
}
impl Gene {
    pub(crate) fn is_reversed(&self) -> bool {
        self.m_reversed
    }

    pub(crate) fn new(m_name: String, m_chr: String, m_start: i32, m_end: i32) -> Self {
        Self {
            m_name,
            m_chr,
            m_start,
            m_end,
            ..Default::default()
        }
    }

    pub(crate) fn valid(&self) -> bool {
        self.m_name != "invalid" && self.m_start != 0 && self.m_end != 0
    }

    pub(crate) fn parse(line_str: &str) -> Result<Self, Error> {
        let splitted = line_str.split(",").collect::<Vec<&str>>();
        log::debug!("splitted={:?}", splitted);

        if splitted.len() < 2 {
            return Ok(Gene::default());
        }

        let name = splitted
            .get(0)
            .unwrap()
            .chars()
            .skip(1)
            .take(splitted.get(0).unwrap().chars().count() - 1)
            .collect::<String>()
            .trim()
            .to_string();

        let chr_pos = splitted.get(1).unwrap().split(":").collect::<Vec<&str>>();
        if chr_pos.len() < 2 {
            return Ok(Gene::default());
        }

        let chr = chr_pos.get(0).unwrap().trim().to_string();

        let range = chr_pos.get(1).unwrap().split("-").collect::<Vec<&str>>();
        if range.len() < 2 {
            return Ok(Gene::default());
        }

        let start = range
            .get(0)
            .unwrap()
            .trim()
            .parse::<i32>()
            // .or_exp()?;
            .or_else(|e| Err(format!("{e:?}, {}", range.get(0).unwrap().trim())))?;
        let end = range
            .get(1)
            .unwrap()
            .trim()
            .parse::<i32>()
            // .or_exp()?;
            .or_else(|e| Err(format!("{e:?}, {}", range.get(1).unwrap().trim())))?;

        Ok(Gene::new(name, chr, start, end))
    }

    pub(crate) fn add_exon(&mut self, id: i32, start: i32, end: i32) -> () {
        let exon = Exon { id, start, end };

        self._add_exon(exon);
    }

    fn _add_exon(&mut self, exon: Exon) {
        let m_exons = &mut self.m_exons;

        m_exons.push(exon);
        if m_exons.len() > 1 {
            if m_exons.get(0).unwrap().start > m_exons.get(1).unwrap().start {
                self.m_reversed = true;
            }
        }
    }

    pub(crate) fn print(&self) -> () {
        println!(
            "{},{}:{}-{}",
            self.m_name, self.m_chr, self.m_start, self.m_end
        );
        println!("{}", {
            if self.m_reversed {
                " reversed"
            } else {
                " forward"
            }
        });

        for i in (0..self.m_exons.len()) {
            println!(
                "{},{},{}",
                self.m_exons.get(i).unwrap().id,
                self.m_exons.get(i).unwrap().start,
                self.m_exons.get(i).unwrap().end
            );
        }
    }

    pub(crate) fn pos2str(&self, pos: i32) -> Result<String, Error> {
        let pp = pos.abs() + self.m_start;

        let mut ss = format!("{}:", self.m_name);

        for i in (0..self.m_exons.len()) {
            if pp >= self.m_exons.get(i).unwrap().start && pp <= self.m_exons.get(i).unwrap().end {
                write!(&mut ss, "exon:{}|", self.m_exons.get(i).unwrap().id)?;
                break;
            }

            if i > 0 {
                if self.m_reversed {
                    if self.m_exons.get(i).unwrap().end < pp
                        && pp < self.m_exons.get(i - 1).unwrap().start
                    {
                        write!(&mut ss, "intron:{}|", self.m_exons.get(i).unwrap().id - 1)?;
                        break;
                    }
                } else {
                    if self.m_exons.get(i - 1).unwrap().end < pp
                        && pp < self.m_exons.get(i).unwrap().start
                    {
                        write!(&mut ss, "intron:{}|", self.m_exons.get(i).unwrap().id - 1)?;
                        break;
                    }
                }
            }
        }

        if pos >= 0 {
            ss.push_str("+");
        } else {
            ss.push_str("-");
        }

        write!(&mut ss, "{}:{}", self.m_chr, pp,)?;

        Ok(ss)
    }

    pub(crate) fn get_exon_intron(
        &self,
        pos: i32,
        is_exon: &mut bool,
        number: &mut i32,
    ) -> () {
        let pp = pos.abs() + self.m_start;

        let mut prev_exon:Option<&Exon> = self.m_exons.first();
        for (i, exon) in self.m_exons.iter().enumerate() {
            if pp>= exon.start && pp <= exon.end {
                *is_exon = true;
                *number = exon.id;
                break;
            }

            if i > 0 {
                if self.m_reversed {
                    if exon.end < pp && pp < prev_exon.unwrap().start {
                        *is_exon=false;
                        *number = exon.id-1;
                        break;
                    }
                } else {
                    if prev_exon.unwrap().end < pp && pp < exon.start {
                        *is_exon = false;
                        *number = exon.id-1;
                        break;
                    }
                }
            }
        }
    }

    pub(crate) fn gene_pos_2_chr_pos(&self, genepos: i32) -> i32 {
        let mut chrpos = genepos.abs() + self.m_start;
        if genepos < 0 {
            chrpos *= -1;
        }

        chrpos
    }
    
}

impl Default for Gene {
    fn default() -> Self {
        Self {
            m_name: "invalid".to_string(),
            m_chr: "invalid".to_string(),
            m_start: 0,
            m_end: 0,
            m_exons: Default::default(),
            m_reversed: false,
        }
    }
}
