use std::ops::DerefMut;
use std::{error::Error, fs::File, io::BufWriter};
use std::io::Write;

use chrono::Local;

use crate::aux::global_settings::global_settings;
use crate::core::html_reporter::FUSIONSCAN_VER;
use crate::genefuse::COMMAND;

use super::{fusion_mapper::FusionMapper, fusion_result::FusionResult};

pub(crate) struct JsonReporter<F> {
    m_filename: String,
    m_fusion_mapper: F,
    m_file: BufWriter<File>,
    // m_fusion_results: Vec<FusionResult>,
}

impl<F> JsonReporter<F>
where
    F: DerefMut<Target = FusionMapper>
{
    pub(crate) fn new(file_name:String, mapper:F) -> Result<Self, Box<dyn Error>> {
        let f = BufWriter::new(File::create(&file_name)?);
        Ok(Self {
            m_filename: file_name,
            m_fusion_mapper: mapper,
            m_file: f,
            // m_fusion_results: ,
        })
    }

    pub(crate) fn run(&mut self) -> Result<(), Box<dyn Error>>{
        let f = &mut self.m_file;

        writeln!(f, "{{", )?;
        writeln!(f, "\t\"command\":\"{}\",", COMMAND.get().unwrap())?;
        writeln!(f, "\t\"version\":\"{}\",", FUSIONSCAN_VER)?;
        writeln!(f, "\t\"time\":\"{}\",", Local::now())?;
        write!(f, "\t\"fusions\":{{")?;

        let mut is_first_mut = true;

        for (i, fusion) in self.m_fusion_mapper.m_fusion_results.iter().enumerate() {
            let matches = fusion.m_matches.as_slice();
            if !global_settings().output_deletions && fusion.is_deletion() {
                continue;
            }
            if fusion.is_left_protein_forward() != fusion.is_right_protein_forward() {
                if !global_settings().output_untranslated {
                    continue;
                }
            }

            if is_first_mut {
                writeln!(f, "")?;
                is_first_mut = false;
            } else {
                write!(f, ",\n")?;
            }

            writeln!(f, "\t\t\"{}\":{{", fusion.m_title)?;
                writeln!(f, "\t\t\t\"left\":{{",)?;
                    writeln!(f, "\t\t\t\t\"gene_name\":\"{}\",", fusion.m_left_gene.m_name)?;
                    writeln!(f, "\t\t\t\t\"gene_chr\":\"{}\",", fusion.m_left_gene.m_chr)?;
                    writeln!(f, "\t\t\t\t\"position\":{},", fusion.m_left_gene.gene_pos_2_chr_pos(fusion.m_left_gp.position))?; //fusion.mLeftGene.genePos2ChrPos(fusion.mLeftGP.position)
                    writeln!(f, "\t\t\t\t\"reference\":\"{}\",", fusion.m_left_ref)?;
                    writeln!(f, "\t\t\t\t\"ref_ext\":\"{}\",", fusion.m_left_ref_ext)?;
                    writeln!(f, "\t\t\t\t\"pos_str\":\"{}\",", fusion.m_left_pos)?;
                    writeln!(f, "\t\t\t\t\"exon_or_intron\":\"{}\",", {if fusion.m_left_is_exon {"exon"} else {"intron"}})?; //(fusion.mLeftIsExon?"exon":"intron")
                    writeln!(f, "\t\t\t\t\"exon_or_intron_id\":{},",fusion.m_left_exon_or_intron_id)?;
                    writeln!(f, "\t\t\t\t\"strand\":\"{}\"", {if fusion.is_left_protein_forward(){"forward"}else{"reversed"}})?;
                writeln!(f, "\t\t\t}}, ",)?;
                writeln!(f, "\t\t\t\"right\":{{",)?;
                    writeln!(f, "\t\t\t\t\"gene_name\":\"{}\",",fusion.m_right_gene.m_name)?;
                    writeln!(f, "\t\t\t\t\"gene_chr\":\"{}\",",fusion.m_right_gene.m_chr)?;
                    writeln!(f, "\t\t\t\t\"position\":{},", fusion.m_right_gene.gene_pos_2_chr_pos(fusion.m_right_gp.position))?;
                    writeln!(f, "\t\t\t\t\"reference\":\"{}\",", fusion.m_right_ref)?;
                    writeln!(f, "\t\t\t\t\"ref_ext\":\"{}\",", fusion.m_right_ref_ext)?;
                    writeln!(f, "\t\t\t\t\"pos_str\":\"{}\",", fusion.m_right_pos)?;
                    writeln!(f, "\t\t\t\t\"exon_or_intron\":\"{}\",", {if fusion.m_right_is_exon {"exon"} else {"intron"}})?; 
                    writeln!(f, "\t\t\t\t\"exon_or_intron_id\":{},", fusion.m_right_exon_or_intron_id)?;
                    writeln!(f, "\t\t\t\t\"strand\":\"{}\"", {if fusion.is_right_protein_forward(){"forward"}else{"reversed"}})?;
                writeln!(f, "\t\t\t}}, ",)?;

            writeln!(f, "\t\t\t\"unique\":{},", fusion.m_unique)?;
            writeln!(f, "\t\t\t\"reads\":[",)?;

            for (m, me) in matches.iter().enumerate() {
                writeln!(f, "\t\t\t\t{{",)?;
                writeln!(f, "\t\t\t\t\t\"break\":{},", me.m_read_break)?;
                writeln!(f, "\t\t\t\t\t\"strand\":\"{}\",", {if me.m_reversed {"reversed"} else {"forward"}})?;
                me.print_read_to_json(f, "\t\t\t\t\t")?;
                write!(f, "\t\t\t\t}}")?;

                if m != matches.len()-1 {
                    write!(f, ",")?;
                }

                write!(f, "\n")?;
            }

            writeln!(f, "\t\t\t]",)?;
            write!(f, "\t\t}}")?;

        }

        writeln!(f, "\n\t}}\n}}\n")?;

        Ok(())
    }




}