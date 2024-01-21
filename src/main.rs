use core::fusion::Fusion;

use argparse::set_configs;
use genefuse::genefuse;

mod argparse;
mod genefuse;
mod aux;
mod utils;
mod core;


#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

fn main() {
    let config = set_configs();

    genefuse(config)
}
