use std::{
    borrow::BorrowMut,
    env::{self},
    process::exit,
    time::Instant,
};

use crate::{
    argparse::RunConfig,
    aux::global_settings::{global_settings, global_settings_w},
    core::fusion_scan::FusionScan,
    utils::check_file_valid,
};

pub(crate) fn genefuse(config: RunConfig) {
    prepare_run(&config);

    log::debug!("start with {} threads", config.thread_num,);
    let timer = Instant::now();

    let fs = FusionScan::new(
        config.fusion_file,
        config.ref_file,
        config.r1_file,
        config.r2_file,
        config.html,
        config.json,
        config.thread_num,
    );

    

    todo!()
}

fn prepare_run(config: &RunConfig) {
    {
        let mut global_settings = global_settings_w();
        global_settings.set_unique_requirement(config.unique);
        global_settings.set_deletion_threshold(config.deletion);
        global_settings.set_output_deletions(config.output_deletion);
        global_settings.set_output_untranslated(config.output_untranslated);
    }

    if config.ref_file.ends_with(".gz") || config.ref_file.ends_with(".gz") {
        println!(
            "reference fasta file should not be compressed.\n\
            please unzip {} and try again.",
            config.ref_file
        );
        exit(-1);
    }

    let command = env::args()
        .into_iter()
        .reduce(|mut a, b| {
            a.push_str(" ");
            a.push_str(&b);
            a
        })
        .unwrap();

    check_file_valid(&config.ref_file);
    check_file_valid(&config.r1_file);

    if config.r2_file != "" {
        check_file_valid(&config.r2_file);
    }

    if config.fusion_file != "" {
        check_file_valid(&config.fusion_file);
    }

    println!("\n# {}\n", command);
}
