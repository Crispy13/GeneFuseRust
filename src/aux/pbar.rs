use std::fmt::Write;

use indicatif::{ProgressBar, ProgressDrawTarget, ProgressState, ProgressStyle};

use crate::core::fusion_scan::MULTI_CSV_MODE;

pub fn prepare_pbar(len: u64) -> ProgressBar {
    // will not use progress bar when in multi-csv-mode.
    if *MULTI_CSV_MODE.get().unwrap() {
        return ProgressBar::hidden();
    }

    _prepare_pbar(len)
    
}

/// Make a progress bar ignoring suppresing conditions. (e.g. multi-csv-mode)
pub(crate) fn prepare_pbar_force(len: u64) -> ProgressBar {
    _prepare_pbar(len)
}

fn _prepare_pbar(len: u64) -> ProgressBar {
    let pb = ProgressBar::new(len);
    pb.set_draw_target(ProgressDrawTarget::stderr_with_hz(8));

    let template = match len {
        1.. => {
            "{spinner:.green} [{elapsed_precise}] {msg} [{bar:.cyan/blue}] {pos}/{len} ({eta}, {per_sec})"
            // "{spinner:.green} [{elapsed_precise}] {msg} [{bar}] {pos}/{len} ({eta}, {per_sec})"
        },
        0 => {
            "{spinner:.green} [{elapsed_precise}] {msg} [ ? ] {pos} ({per_sec})"
        }
    };
    
    pb.set_style(
        ProgressStyle::with_template(
            template
        )
        .unwrap()
        .with_key("eta", |state: &ProgressState, w: &mut dyn Write| {
            write!(w, "{:.1}s", state.eta().as_secs_f64()).unwrap()
        })
        // .progress_chars("#>-"),
    );

    pb
}

pub trait PBSummary {
    fn finish_with_summary(&self);

    fn finish_with_summary_force(&self);
}

impl PBSummary for ProgressBar {
    fn finish_with_summary(&self) {
        if self.is_hidden() {
            let elapsed_secs = self.elapsed().as_secs_f64();
            let pb_pos = self.position();
            eprintln!(
                "[{elapsed_precise}] {pos} ({per_sec:.2}/s)",
                elapsed_precise = get_hms(elapsed_secs),
                pos = pb_pos,
                per_sec = pb_pos as f64 / elapsed_secs,
            )
        }

        self.finish();
    }

    fn finish_with_summary_force(&self) {
        let elapsed_secs = self.elapsed().as_secs_f64();
        let pb_pos = self.position();
        eprintln!(
            "[{elapsed_precise}] {pos} ({per_sec:.2}/s)",
            elapsed_precise = get_hms(elapsed_secs),
            pos = pb_pos,
            per_sec = pb_pos as f64 / elapsed_secs,
        );

        self.finish_and_clear();
    }
}

#[inline]
pub fn get_hms(dur_secs:f64) -> String {
    let (hours, rem) = div_mod(dur_secs, 3600_f64);
    let (mins, rem) = div_mod(rem, 60_f64);
    let secs = rem % 60_f64;

    format!("{:0>2}:{:0>2}:{:0>4.1}", hours, mins, secs)
}

fn div_mod(n: f64, d: f64) -> (f64, f64) {
    ((n / d).trunc(), n % d)
}