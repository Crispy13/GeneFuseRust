use log::LevelFilter;
use log4rs::append::console::ConsoleAppender;
use log4rs::append::file::FileAppender;
use log4rs::encode::pattern::PatternEncoder;
use log4rs::config::{Appender, Config, Logger, Root};

pub(crate) fn set_logger() -> log4rs::Handle {
    let stderr = ConsoleAppender::builder().target(log4rs::append::console::Target::Stderr).build();

    let pattern_encoder = Box::new(PatternEncoder::new("[{d}] {l}> - {m}{n}"));


    // let requests = FileAppender::builder()
    //     .encoder(Box::new(PatternEncoder::new("[{d}] {l}> - {m}{n}")))
    //     .build("log/requests.log")
    //     .unwrap();

    let config = Config::builder()
        .appender(Appender::builder().build("stderr", Box::new(stderr)))
        // .logger(Logger::builder().build("app::backend::db", LevelFilter::Info))
        // .logger(Logger::builder()
        //     .appender("requests")
        //     .additive(false)
        //     .build("app::requests", LevelFilter::Info))
        .build(Root::builder().appender("stderr").build(LevelFilter::Debug))
        .unwrap();

    log4rs::init_config(config).unwrap()
}


