use clap::{arg, command, value_parser, ArgAction, ArgMatches};

pub(crate) fn parse_args() -> ArgMatches {
    let command = command!() // requires `cargo` feature
        .arg(
            arg!(
                -'1' --read1 <read1> "read1 file name"
            )
            .required(true)
            .value_parser(value_parser!(String)),
        )
        .arg(
            arg!(
                -'2' --read2 <read2> "read2 file name"
            )
            .required(false)
            .value_parser(value_parser!(String))
            .default_value("")
        )
        .arg(
            arg!(
                -f --fusion <fusion> "fusion file name, in CSV format"
            )
            .required(true)
            .value_parser(value_parser!(String)),
        )
        .arg(
            arg!(
                -r --ref <ref> "reference fasta file name"
            )
            .required(true)
            .value_parser(value_parser!(String)),
        )
        .arg(
            arg!(
                -u --unique <unique> "specify the least supporting read number is required to report a fusion, default is 2"
            )
            .required(false)
            .value_parser(value_parser!(usize))
            .default_value("2")
        )
        .arg(
            arg!(
                -h --html <html> "file name to store HTML report, default is genefuse.html"
            )
            .required(false)
            .value_parser(value_parser!(String))
            .default_value("genefuse.html")
        )
        .arg(
            arg!(
                -j --json <json> "file name to store JSON report, default is genefuse.json"
            )
            .required(false)
            .value_parser(value_parser!(String))
            .default_value("genefuse.json")
        )
        .arg(
            arg!(
                -t --thread <thread> "worker thread number, default is 4"
            )
            .required(false)
            .value_parser(value_parser!(usize))
            .default_value("4")
        )
        .arg(
            arg!(
                -d --deletion <deletion> "specify the least deletion length of a intra-gene deletion to report, default is 50"
            )
            .required(false)
            .value_parser(value_parser!(usize))
            .default_value("50")
        )
        .arg(
            arg!(
                -D --output_deletions <output_deletions> "long deletions are not output by default, enable this option to output them"
            )
            .required(false)
            .action(ArgAction::SetTrue)
        )
        .arg(
            arg!(
                -U --output_untranslated_fusions <output_untranslated_fusions> 
            )
            .help("the fusions that cannot be transcribed or translated are not output by default, \
            enable this option to output them")
            .required(false)
            .action(ArgAction::SetTrue)
        );

    command.get_matches()
}

pub(crate) struct RunConfig {
    pub(crate) r1_file: String,
    pub(crate) r2_file: String,
    pub(crate) fusion_file: String,
    pub(crate) html: String,
    pub(crate) json: String,
    pub(crate) ref_file: String,
    pub(crate) thread_num: usize,
    pub(crate) unique: usize,
    pub(crate) deletion: usize,
    pub(crate) output_deletion: bool,
    pub(crate) output_untranslated: bool,
}

impl RunConfig {
    fn from_args(mut args: ArgMatches) -> RunConfig {
        Self {
            r1_file: args.remove_one::<String>("read1").unwrap(),
            r2_file: args.remove_one::<String>("read2").unwrap(),
            fusion_file: args.remove_one::<String>("fusion").unwrap(),
            html: args.remove_one::<String>("html").unwrap(),
            json: args.remove_one::<String>("json").unwrap(),
            ref_file: args.remove_one::<String>("ref").unwrap(),
            thread_num: args.remove_one::<usize>("thread").unwrap(),
            unique: args.remove_one::<usize>("unique").unwrap(),
            deletion: args.remove_one::<usize>("deletion").unwrap(),
            output_deletion: args.remove_one::<bool>("output_deletions").unwrap(),
            output_untranslated: args.remove_one::<bool>("output_untranslated_fusions").unwrap(),
        }
    
    }
}

pub(crate) fn set_configs() -> RunConfig {

    RunConfig::from_args(parse_args())
}
