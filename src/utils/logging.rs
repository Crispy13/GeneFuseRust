use log::LevelFilter;
use log4rs::append::console::ConsoleAppender;
use log4rs::append::file::FileAppender;
use log4rs::config::{Appender, Config, Logger, Root};
use log4rs::encode::pattern::PatternEncoder;

pub(crate) fn init_logger() -> log4rs::Handle {
    let pattern_encoder = Box::new(PatternEncoder::new("[{d}] {t} {l}> - {m}{n}"));

    let stderr = ConsoleAppender::builder()
        .target(log4rs::append::console::Target::Stderr)
        .encoder(Box::new(PatternEncoder::new(
            "[{d(%Y-%m-%d %H:%M:%S)}] {T} {t} {l}>> {m}{n}",
        )))
        .build();

    // let file_appender = FileAppender::builder()
    //     .encoder(pattern_encoder)
    //     .append(false)
    //     .build("debug.log")
    //     .unwrap();

    let config = Config::builder()
        // .appender(Appender::builder().build("file", Box::new(file_appender)))
        .appender(Appender::builder().build("stderr", Box::new(stderr)))
        // .logger(Logger::builder().build("app::backend::db", LevelFilter::Info))
        // .logger(Logger::builder()
        //     .appender("requests")
        //     .additive(false)
        //     .build("app::requests", LevelFilter::Info))
        .build(
            Root::builder()
                // .appender("file")
                .appender("stderr")
                .build(LevelFilter::Info),
        )
        .unwrap();

    log4rs::init_config(config).unwrap()
}
