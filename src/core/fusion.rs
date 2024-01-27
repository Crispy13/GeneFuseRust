use std::{
    default,
    error::Error,
    fs::File,
    io::{BufRead, BufReader, Read},
};

use super::gene::Gene;

// use anyhow::Result;
#[derive(Debug, Clone)]
pub(crate) struct Fusion {
    pub(crate) m_gene: Gene,
}

impl Fusion {
    pub(crate) fn new(gene: Gene) -> Fusion {
        Self { m_gene: gene }
    }

    pub(crate) fn parse_csv(filename: &str) -> Result<Vec<Self>, Box<dyn Error>> {
        const max_line: usize = 4096;

        let mut file = BufReader::new(File::open(filename)?).take(max_line as u64);
        let mut line = Vec::<char>::with_capacity(max_line);

        let mut fusions: Vec<Self> = Vec::new();
        let mut working_gene = Gene::default();

        let mut line_s = String::new();
        while let Ok(true) = {
            line.clear();
            line_s.clear();
            file.read_line(&mut line_s).and_then(|rl| Ok(rl > 0))
        } {
            // trim \n, \r or \r\n in the tail
            line.extend(line_s.chars());
            let readed = line.len();

            if readed >= 2 {
                if line.get(readed - 1).unwrap().eq(&'\n')
                    || line.get(readed - 1).unwrap().eq(&'\r')
                {
                    *line.get_mut(readed - 1).unwrap() = '\0';
                    if line.get(readed - 2).unwrap().eq(&'\r') {
                        *line.get_mut(readed - 2).unwrap() = '\0';
                    }
                }
            }

            let line_str = line_s.trim();

            let splitted = line_str.split(",").collect::<Vec<&str>>();
            // wrong line
            if splitted.len() < 2 {
                continue;
            }
            // comment line
            if splitted.first().unwrap().starts_with("#") {
                continue;
            }
            // gene line
            if splitted.first().unwrap().starts_with(">") {
                if working_gene.valid() {
                    fusions.push(Fusion::new(working_gene));
                }

                working_gene = Gene::parse(line_str)?;
                continue;
            }

            // position line require id, start, position
            if splitted.len() < 3 {
                continue;
            }

            let id = splitted.get(0).unwrap().trim().parse::<i32>()?;
            let start = splitted.get(1).unwrap().trim().parse::<i32>()?;
            let end = splitted.get(2).unwrap().trim().parse::<i32>()?;
            working_gene.add_exon(id, start, end);
        }
        // last one
        if working_gene.valid() {
            fusions.push(Fusion::new(working_gene));
        }

        Ok(fusions)
    }

    fn print(&self) -> () {
        self.m_gene.print()
    }

    fn print_html(&self, file: File) -> () {
        todo!()
    }

    pub(crate) fn is_reversed(&self) -> bool {
        self.m_gene.is_reversed()
    }

    fn pos2str(&self, pos: i32) -> Result<String, Box<dyn Error>> {
        self.m_gene.pos2str(pos)
    }
}

#[cfg(test)]
mod test {
    use super::Fusion;
    use anyhow::Result;

    #[test]
    fn test1() {
        let fusions = Fusion::parse_csv("testdata/fusions.csv").unwrap();

        let _inner = || {
            for i in (0..fusions.len()) {
                let f = fusions.get(i).unwrap();
                if f.m_gene.m_name == "ALK" {
                    // exon 20
                    if f.pos2str(-30582).unwrap() != "ALK:exon:20|-chr2:29446222" {
                        return f.pos2str(-30582).unwrap();
                    }
                    // intron 19
                    if f.pos2str(31060).unwrap() != "ALK:intron:19|+chr2:29446700" {
                        return f.pos2str(31060).unwrap();
                    }
                }

                if f.m_gene.m_name == "EML4" {
                    // exon6
                    if f.pos2str(95365).unwrap() != "EML4:exon:6|+chr2:42491855" {
                        return f.pos2str(95365).unwrap();
                    }
                    // intron 5
                    if (f.pos2str(95346).unwrap() != "EML4:intron:5|+chr2:42491836") {
                        return f.pos2str(95346).unwrap();
                    }
                }
            }

            return "true".to_string();
        };

        println!("{}", _inner());
    }
}
