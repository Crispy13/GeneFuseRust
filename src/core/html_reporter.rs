use std::{
    error, fs::File, io::{BufWriter, Write}, ops::DerefMut
};

use chrono::Local;

use crate::{aux::global_settings::global_settings, genefuse::COMMAND};

use super::{fusion_mapper::FusionMapper, fusion_result::FusionResult, fusion_scan::Error};

pub(crate) const FUSIONSCAN_VER: &str = env!("CARGO_PKG_VERSION");

pub(crate) struct HtmlReporter<'f, 's> {
    m_filename: String,
    m_fusion_mapper: &'f mut FusionMapper<'s>,
    m_file: BufWriter<File>,
    // m_fusion_results: Vec<FusionResult>,
}

impl<'f, 's> HtmlReporter<'f, 's>
{
    pub(crate) fn new(
        filename: String,
        mapper: &'f mut FusionMapper<'s>,
    ) -> Result<Self, Error> {
        let m_file = BufWriter::new(File::create(&filename)?);
        Ok(Self {
            m_filename: filename,
            m_fusion_mapper: mapper,
            m_file: m_file,
            // m_fusion_results: mapper.m_fusion_results,
        })
    }

    fn m_fusion_results(&self) -> &[FusionResult] {
        &self.m_fusion_mapper.m_fusion_results
    }

    pub(crate) fn run(&mut self) -> Result<(), Error> {
        log::debug!("printing header...");
        self.print_header()?;
        log::debug!("printing helper...");
        self.print_helper()?;
        log::debug!("printing fusions...");
        self.print_fusions()?;
        log::debug!("printing footer...");
        self.print_footer()?;

        Ok(())
    }

    fn print_header(&mut self) -> Result<(), Error> {
        let f = &mut self.m_file;

        write!(
            f,
            "<html><head><meta http-equiv=\"content-type\" content=\"text/html;charset=utf-8\" />"
        )?;
        write!(
            f,
            "<title>GeneFuse {}, at {}</title>",
            FUSIONSCAN_VER,
            Local::now()
        )?;

        self.print_js();
        self.print_css();

        let f = &mut self.m_file;

        write!(f, "</head>")?;
        write!(f, "<body><div id='container'>")?;
        write!(
            f,
            "<div class='software'> \
                <a href='https://github.com/OpenGene/GeneFuse' style='text-decoration:none;' \
                target='_blank'>GeneFuse</a> <font size='-1'>{}</font></div>",
            FUSIONSCAN_VER
        )?;

        Ok(())
    }

    fn print_css(&mut self) -> Result<(), Error> {
        let f = &mut self.m_file;

        write!(f, "<style type=\"text/css\">")?;
        write!(
            f,
            "td {{border:1px solid #dddddd;padding-left:2px;padding-right:2px;font-size:10px;}}"
        )?;
        write!(
            f,
            "table {{border:1px solid #999999;padding:2x;border-collapse:collapse;}}"
        )?;
        write!(f, "img {{padding:30px;}}")?;
        write!(f, ".alignleft {{text-align:left;}}")?;
        write!(f, ".alignright {{text-align:right;}}")?;
        write!(
            f,
            ".software {{font-weight:bold;font-size:24px;padding:5px;}}"
        )?;
        write!(
            f,
            ".header {{color:#ffffff;padding:1px;height:20px;background:#000000;}}"
        )?;
        write!(
            f,
            ".figuretitle {{color:#996657;font-size:20px;padding:50px;}}"
        )?;
        write!(f, "#container {{text-align:center;padding:1px;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}}")?;
        write!(
            f,
            "#menu {{padding-top:10px;padding-bottom:10px;text-align:left;}}"
        )?;
        write!(f, "#menu a {{color:#0366d6; font-size:18px;font-weight:600;line-height:28px;text-decoration:none;font-family:-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol'}}")?;
        write!(f, "a:visited {{color: #999999}}")?;
        write!(
            f,
            ".menu_item {{text-align:left;padding-top:5px;font-size:18px;}}"
        )?;
        write!(f, ".highlight {{text-align:left;padding-top:30px;padding-bottom:30px;font-size:20px;line-height:35px;}}")?;
        write!(f, ".fusion_head {{text-align:left;color:#0092FF;font-family:Arial;padding-top:20px;padding-bottom:5px;}}")?;
        write!(f, ".fusion_block {{}}")?;
        write!(f, ".match_brief {{font-size:8px}}")?;
        write!(f, ".fusion_point {{color:#FFCCAA}}")?;
        write!(
            f,
            "#helper {{text-align:left;border:1px dotted #fafafa;color:#777777;font-size:12px;}}"
        )?;
        write!(f, "#footer {{text-align:left;padding-left:10px;padding-top:20px;color:#777777;font-size:10px;}}")?;
        write!(
            f,
            ".exon_left{{background:blue;color:white;border:0px;padding:0px;font-size:8px;}}"
        )?;
        write!(
            f,
            ".exon_right{{background:red;color:white;0px;padding:0px;font-size:8px;}}"
        )?;
        write!(
            f,
            ".intron_left{{color:blue;0px;padding:0px;font-size:8px;}}"
        )?;
        write!(
            f,
            ".intron_right{{color:red;0px;padding:0px;font-size:8px;}}"
        )?;
        write!(f, ".protein_table{{text-align:center;font-size:8px;}}")?;
        write!(
            f,
            ".tips{{font-size:10px;padding:5px;color:#666666;text-align:left;}}"
        )?;

        write!(f, "</style>")?;

        Ok(())
    }

    fn print_js(&mut self) -> Result<(), Error> {
        let f = &mut self.m_file;

        write!(f, "<script type=\"text/javascript\">\n")?;
        write!(
            f,
            "function toggle(targetid){{ \n\
                        if (document.getElementById){{ \n\
                            target=document.getElementById(targetid); \n\
                                if (target.style.display=='table-row'){{ \n\
                                    target.style.display='none'; \n\
                                }} else {{ \n\
                                    target.style.display='table-row'; \n\
                                }} \n\
                        }} \n\
                    }}"
        )?;
        write!(
            f,
            "function toggle_target_list(targetid){{ \n\
                        if (document.getElementById){{ \n\
                            target=document.getElementById(targetid); \n\
                                if (target.style.display=='block'){{ \n\
                                    target.style.display='none'; \n\
                                    document.getElementById('target_view_btn').value='view';\n\
                                }} else {{ \n\
                                    document.getElementById('target_view_btn').value='hide';\n\
                                    target.style.display='block'; \n\
                                }} \n\
                        }} \n\
                    }}"
        )?;
        write!(f, "</script>")?;
        Ok(())
    }

    fn print_footer(&mut self) -> Result<(), Error> {
        let f = &mut self.m_file;
        write!(f, "<div id='footer'> ")?;
        write!(f, "<p>{}</p>", COMMAND.get().unwrap())?;

        self.print_scan_targets()?;

        let f = &mut self.m_file;
        write!(f, "GeneFuse {}, at {} </div>", FUSIONSCAN_VER, Local::now())?;
        write!(f, "</div></body></html>")?;

        Ok(())
    }

    fn print_helper(&mut self) -> Result<(), Error> {
        let f = &mut self.m_file;
        write!(f, "<div id='helper'><p>Helpful tips:</p><ul>",)?;
        write!(f, "<li> Base color indicates quality: <font color='#78C6B9'>extremely high (Q40+)</font>, <font color='#33BBE2'>high (Q30~Q39) </font>, <font color='#666666'>moderate (Q20~Q29)</font>, <font color='#E99E5B'>low (Q15~Q19)</font>, <font color='#FF0000'>extremely low (0~Q14).</font> </li>",)?;
        write!(
            f,
            "<li> Move mouse over the base, it will show the quality value</li>",
        )?;
        write!(
            f,
            "<li> Click on any row, the original read/pair will be displayed</li>",
        )?;
        write!(f, "<li> For pair-end sequencing, GeneFuse tries to merge each pair, with overlapped assigned higher qualities </li>",)?;
        write!(f, "</ul><p>Columns:</p><ul>",)?;
        write!(f, "<li> col1: is fusion mapped with original read? → means original read, ← means reverse complement</li>",)?;
        write!(f, "<li> col2: edit distance (ed) between read and reference sequence (left_part_ed | right_part_ed)</li>",)?;
        write!(f, "<li> col3: read's left part after fusion break</li>",)?;
        write!(f, "<li> col4: read's right part after fusion break</li>",)?;
        write!(f, "</ul></div>",)?;
        Ok(())
    }

    fn print_fusions(&mut self) -> Result<(), Error> {
        // calculate the found fusion
        let found = self.m_fusion_results().len();

        let f = &mut self.m_file;

        // print menu
        write!(f, "<div id='menu'><p>Found {} fusion", found)?;
        if found > 1 {
            write!(f, "s",)?;
        }

        write!(f, ":</p><ul>")?;

        let mut id = 0;
        for fr in self.m_fusion_mapper.m_fusion_results.iter() {
            id += 1;
            write!(
                f,
                "<li class='menu_item'><a href='#fusion_id_{}'> {}, {}</a></li>",
                id, id, fr.m_title
            )?;
            
        }

        write!(f, "</ul></div>")?;

        let mut id = 0;
        log::debug!("print fusion for loop");
        log::debug!("found={}, self.m_fusion_mapper.m_fusion_results.len()={}", found, self.m_fusion_mapper.m_fusion_results.len());

        for mut fusion in self.m_fusion_mapper.m_fusion_results.iter_mut() {
            if !global_settings().output_deletions && fusion.is_deletion() {
                continue;
            }
            if fusion.is_left_protein_forward() != fusion.is_right_protein_forward() {
                if !global_settings().output_untranslated {
                    continue;
                }
            }

            id += 1;
            Self::print_fusion(id, &mut fusion, f)?;
        }

        Ok(())
    }

    fn print_fusion(
        id: i32,
        fusion: &mut FusionResult,
        f: &mut BufWriter<File>,
    ) -> Result<(), Error> {
        // let f = &mut self.m_file;

        

        write!(f, "<div class='fusion_block'>")?;
        write!(f, "<div class='fusion_head'><a name='fusion_id_{}'>", id)?;
        write!(f, "{}, {}", id, fusion.m_title)?;
        write!(f, "</a></div>")?;

        write!(f, "<div class='tips'>Inferred protein")?;

        if fusion.is_left_protein_forward() != fusion.is_right_protein_forward() {
            write!(
                f,
                " (transcription direction conflicts, this fusion may be not transcribed) "
            )?;
        }

        write!(f, ":</div>")?;

        fusion.print_fusion_protein_html(f)?;

        write!(f, "<div class='tips'>Supporting reads:</div>")?;
        write!(f, "<table>")?;
        write!(f, "<tr class='header'>")?;
        write!(
            f,
            "<td class='alignright' colspan='3'>{} = <font color='yellow'>↓</font></td>",
            fusion.m_left_pos
        )?;
        write!(
            f,
            "<td class='alignleft'><font color='yellow'>↓</font> = {}</td>",
            fusion.m_right_pos
        )?;
        write!(f, "</tr>")?;
        write!(f, "<tr class='header'>")?;
        write!(
            f,
            "<td class='alignright' colspan='3'><a title='{}___{}'>{}</a></td>",
            fusion.m_left_ref, fusion.m_left_ref_ext, fusion.m_left_ref
        )?;
        write!(
            f,
            "<td class='alignleft'><a title='{}___{}'>{}</a></td>",
            fusion.m_right_ref_ext, fusion.m_right_ref, fusion.m_right_ref
        )?;
        write!(f, "</tr>")?;
        
        let read_matches = fusion.m_matches.as_slice();
        log::debug!("print_fusion(): read_matches_len={}", read_matches.len());
        for (m, me) in read_matches.iter().enumerate() {
            let rowid = id*100000 + m as i32;
            write!(f, "<tr onclick='toggle({});'>", rowid)?;
            write!(f, "<td>")?;
            write!(f, "<a title='{}'>", me.m_read.m_name)?;
            // for display alignment
            if (m+1) < 10 {
                write!(f, "0")?;
            }
            if (m+1) < 100 {
                write!(f, "0")?;
            }
            if (m+1) < 1000 {
                write!(f, "0")?;
            }
            write!(f, "{}", m+1)?;
            me.print_html_td(f)?;
            write!(f, "</tr>")?;

            // print a hidden row containing the full read
            write!(f, "<tr id='{}' style='display:none;'>", rowid)?;
            write!(f, "<td colspan='6'><xmp>")?;
            me.print_reads_to_file(f)?;
            write!(f, "</xmp></td>")?;
            write!(f, "</tr>")?;
        }
        write!(f, "</table></div>")?;
        
        Ok(())
    }

    fn print_scan_targets(&mut self) -> Result<(), Error> {
        // original cpp code does nothing.
        Ok(())
    }
}
