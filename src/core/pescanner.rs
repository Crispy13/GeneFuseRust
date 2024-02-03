use rayon::{prelude::*, ThreadPoolBuilder};
use std::{
    error::Error,
    io::Write,
    panic::Location,
    process::exit,
    sync::{
        atomic::{AtomicBool, AtomicUsize, Ordering},
        Condvar, Mutex, RwLock,
    },
    thread::sleep,
    time::Duration,
};

use super::{
    common::{PACK_NUM_LIMIT, PACK_SIZE},
    fasta_reader::{self, FastaReader},
    fastq_reader::FastqReaderPair,
    fusion_mapper::FusionMapper,
    read::SequenceReadPair,
    read_match::ReadMatch,
};
use crate::{
    core::{html_reporter::HtmlReporter, json_reporter::JsonReporter},
    utils::open_csv,
};
use crossbeam::queue::ArrayQueue;

pub const DBT: &'static str = "@NB551106:59:HTFV3BGX2:4:22605:18628:10037";

#[derive(Debug)]
struct ReadPairPack {
    data: Vec<SequenceReadPair>,
    count: i32,
}

// struct CPPConditionVariable {}

struct ReadPairRepository {
    pack_buffer: ArrayQueue<ReadPairPack>,
    read_pos: AtomicUsize,
    write_pos: AtomicUsize,
    read_counter: AtomicUsize,
    // mtx: CPPMutex,
    // read_counter_mtx: CPPMutex,
    // repo_not_full: Condvar,
    // repo_not_empty: Condvar,
}

pub(crate) struct PairEndScanner {
    m_fusion_file: String,
    m_ref_file: String,
    m_read1_file: String,
    m_read2_file: String,
    m_html_file: String,
    m_json_file: String,
    m_repo_o: Option<ReadPairRepository>,
    m_produce_finished: AtomicBool,
    m_thread_num: i32,
    m_fusion_mapper_o: Option<FusionMapper>,
}

impl PairEndScanner {
    pub(crate) fn new(
        fusion_file: String,
        ref_file: String,
        read1_file: String,
        read2_file: String,
        html: String,
        json: String,
        thread_num: i32,
    ) -> Self {
        Self {
            m_fusion_file: fusion_file,
            m_ref_file: ref_file,
            m_read1_file: read1_file,
            m_read2_file: read2_file,
            m_html_file: html,
            m_json_file: json,
            m_repo_o: None,
            m_produce_finished: AtomicBool::new(false),
            m_thread_num: thread_num,
            m_fusion_mapper_o: None,
            // repo_not_full: Condvar::new(),
            // repo_not_empty: Condvar::new(),
        }
    }

    // fn m_fusion_mapper(&self) -> &FusionMapper {
    //     self.m_fusion_mapper_o.as_ref().unwrap()
    // }

    // fn m_fusion_mapper_mut(&mut self) -> &mut FusionMapper {
    //     self.m_fusion_mapper_o.as_mut().unwrap()
    // }

    fn m_repo(&self) -> &ReadPairRepository {
        self.m_repo_o.as_ref().unwrap()
    }

    // fn m_repo_mut(&mut self) -> &mut ReadPairRepository {
    //     self.m_repo_o.as_mut().unwrap()
    // }

    fn init_pack_repository(&mut self) {
        let pack_buffer = ArrayQueue::new(PACK_NUM_LIMIT as usize); //Vec::<ReadPairPack>::with_capacity(PACK_NUM_LIMIT as usize);
        let write_pos = AtomicUsize::new(0);
        let read_pos = AtomicUsize::new(0);
        let read_counter = AtomicUsize::new(0);

        if let Some(m_repo) = self.m_repo_o.as_mut() {
            // let mut m_repo = m_repo_lock.unwrap();

            m_repo.pack_buffer = pack_buffer;
            m_repo.write_pos = write_pos;
            m_repo.read_pos = read_pos;
            m_repo.read_counter = read_counter;
        } else {
            let repo = ReadPairRepository {
                pack_buffer,
                read_pos: read_pos,
                write_pos: write_pos,
                read_counter: read_counter,
                // mtx: CPPMutex {},
                // read_counter_mtx: CPPMutex {},
                // repo_not_full: Condvar::new(),
                // repo_not_empty: Condvar::new(),
            };

            self.m_repo_o = Some(repo);
        }
    }

    fn produce_pack(&self, mut pack: ReadPairPack) {
        log::debug!("Entered produce_pack()");
        // lock m_repo.mtx;

        // while (self.m_repo().write_pos + 1) % PACK_NUM_LIMIT as usize == self.m_repo().read_pos {
        //     // self.m_repo().repo_not_full.wait(lock);
        // }

        // let wp = self.m_repo().write_pos;
        let m_repo = self.m_repo_o.as_ref().unwrap();

        loop {
            match m_repo.pack_buffer.push(pack) {
                Ok(_) => break,
                Err(p) => {
                    // means the buffer is full.
                    pack = p;
                }
            }
        }

        // *self.m_repo_mut().pack_buffer.get_mut(wp).unwrap() = pack;
        // m_repo.write_pos.fetch_add(1, Ordering::Relaxed);
        // log::debug!("Add 1 to write_pos.");

        // let _ = m_repo.write_pos.compare_exchange(
        //     PACK_NUM_LIMIT as usize,
        //     0,
        //     Ordering::Relaxed,
        //     Ordering::Relaxed,
        // ); // ignore error because in this method, occured error means just the atomic value is not the same as `current`.

        // mRepo.repoNotEmpty.notify_all();
        // lock.unlock();
    }

    fn producer_task(&self) -> Result<(), Box<dyn Error>> {
        log::debug!("Entered producer_task()");

        let mut slept = 0;
        let mut data = Vec::<SequenceReadPair>::with_capacity(PACK_SIZE as usize);

        let mut reader = FastqReaderPair::from_paths(&self.m_read1_file, &self.m_read2_file)?;

        let mut count = 0;

        let mut read;
        loop {
            read = reader.read();

            if read.is_none() {
                // the last pack
                break;
            }

            data.push(read.unwrap());
            count += 1;

            // a full pack
            if count == PACK_SIZE {
                let pack = ReadPairPack { data, count };
                self.produce_pack(pack);

                //re-initialize data for next pack
                data = Vec::<SequenceReadPair>::with_capacity(PACK_SIZE as usize);

                // reset count to 0
                count = 0;
                // if pack_buffer is full, sleep a while.
                // while self.m_repo_o.as_ref().unwrap().pack_buffer.is_full()
                // {
                //     slept += 1;
                //     sleep(Duration::from_micros(1000));
                // }
            }
        }

        let pack = ReadPairPack { data, count };
        self.produce_pack(pack);

        // lock self.m_repo().read_counter_mtx;
        log::debug!("producing finished.");
        self.m_produce_finished.store(true, Ordering::Relaxed);
        // lock.unlock();

        // producing tasks has been done, from now try to consume the tasks.
        if rayon::current_num_threads() > 1 {
            while !self.m_repo_o.as_ref().unwrap().pack_buffer.is_empty() {
                self.consume_pack()?;
            }
        }
        

        Ok(())
    }

    pub(crate) fn drop_and_get_back_fasta_reader(self) -> FastaReader {
        self.m_fusion_mapper_o.unwrap().m_indexer.m_reference.unwrap()
    }

    pub(crate) fn scan_per_fusion_csv(
        &mut self,
        fasta_reader: FastaReader,
    ) -> Result<bool, Box<dyn Error>> {
        log::debug!("Entered into scan.");
        self.m_fusion_mapper_o = Some(FusionMapper::from_fasta_reader_and_fusion_files(
            fasta_reader,
            &self.m_fusion_file,
        )?);
        log::debug!("Made fusion mapper.");

        self.init_pack_repository();

        // log::debug!("fusion_list={:#?}", self.m_fusion_mapper_o.as_ref().unwrap().fusion_list);

        // for (k,v) in self.m_fusion_mapper_o.as_ref().unwrap().m_indexer.m_reference.as_ref().unwrap().m_all_contigs.iter() {
        //     log::debug!("contig={}, seq_len={}, seq_5_from_start={}, seq_last_5_to_end={}", k, v.len(), v.get(..5).unwrap(), v.get((v.len()-5)..).unwrap())
        // }

        Ok(self._scan()?)
    }

    pub(crate) fn scan(&mut self) -> Result<bool, Box<dyn Error>> {
        log::debug!("Entered into scan.");
        self.m_fusion_mapper_o = Some(FusionMapper::from_ref_and_fusion_files(
            &self.m_ref_file,
            &self.m_fusion_file,
        )?);
        log::debug!("Made fusion mapper.");

        self.init_pack_repository();

        // log::debug!("fusion_list={:#?}", self.m_fusion_mapper_o.as_ref().unwrap().fusion_list);

        // for (k,v) in self.m_fusion_mapper_o.as_ref().unwrap().m_indexer.m_reference.as_ref().unwrap().m_all_contigs.iter() {
        //     log::debug!("contig={}, seq_len={}, seq_5_from_start={}, seq_last_5_to_end={}", k, v.len(), v.get(..5).unwrap(), v.get((v.len()-5)..).unwrap())
        // }

        Ok(self._scan()?)
    }

    fn _scan(&mut self) -> Result<bool, Box<dyn Error>> {
        rayon::scope(|tps| {
            tps.spawn(|tps| self.producer_task().unwrap());

            for t in (0..(rayon::current_num_threads() - 1)) {
                tps.spawn(|tps| {
                    self.consumer_task().unwrap();
                })
            }
        });

        log::debug!("Produced and consumed all the tasks.");
        let m_fusion_mapper = self.m_fusion_mapper_o.as_mut().unwrap();

        // let mut obj_csv = open_csv();
        // writeln!(&mut obj_csv, "i,j,mName").unwrap();
        // m_fusion_mapper
        //     .fusion_matches
        //     .lock()
        //     .unwrap()
        //     .iter()
        //     .enumerate()
        //     .for_each(|(i, v)| {
        //         v.iter().enumerate().for_each(|(j, rm)| {
        //             writeln!(&mut obj_csv, "{},{},{}", i, j, rm.m_read.m_name).unwrap();
        //         })
        //     });

        // drop(obj_csv);

        // exit(0);

        log::debug!("run matches methods...");
        m_fusion_mapper.filter_matches();
        m_fusion_mapper.sort_matches();
        m_fusion_mapper.cluster_matches();

        log::debug!("making html reports...");
        self.html_report().unwrap();
        log::debug!("making json reports...");
        self.json_report().unwrap();

        let m_fusion_mapper = self.m_fusion_mapper_o.as_mut().unwrap();
        m_fusion_mapper.free_matches();

        Ok(true)
    }

    fn consumer_task(&self) -> Result<(), Box<dyn Error>> {
        log::debug!("Entered consumer_task()");
        loop {
            // lock = self.m_repo().read_counter_mtx;
            if self.m_produce_finished.load(Ordering::Relaxed)
                && self.m_repo_o.as_ref().unwrap().pack_buffer.is_empty()
            {
                // lock.unlock();
                break;
            }

            if self.m_produce_finished.load(Ordering::Relaxed) {
                self.consume_pack()?;
                // lock.unlock();
            } else {
                // lock.unlock();
                self.consume_pack()?;
            }
        }

        log::debug!("consuming finished.");
        Ok(())
    }

    fn consume_pack(&self) -> Result<(), Box<dyn Error>> {
        // let data = ReadPairPack {
        //     data: todo!(),
        //     count: todo!(),
        // };

        // lock self.m_repo().mtx;

        // read buffer is empty, just wait here.
        // while self.m_repo().write_pos.load(Ordering::Relaxed) as i32 % PACK_NUM_LIMIT
        //     == self.m_repo().read_pos.load(Ordering::Relaxed) as i32 % PACK_NUM_LIMIT
        // {
        //     if self.m_produce_finished.load(Ordering::Relaxed) {
        //         // lock.unlock();
        //         return Ok(());
        //     }
        //     // self.m_repo().repo_not_empty.wait(lock);
        // }

        // let data = self
        //     .m_repo_o
        //     .as_ref()
        //     .unwrap()
        //     .pack_buffer
        //     .get(self.m_repo().read_pos)
        //     .unwrap();

        let data = loop {
            if let Some(data) = self.m_repo_o.as_ref().unwrap().pack_buffer.pop() {
                break data;
            } else {
                if self.m_produce_finished.load(Ordering::Relaxed) {
                    return Ok(());
                }
            }
        };

        // let data = self.m_repo_o.as_ref().unwrap().pack_buffer.pop().unwrap();

        // self.m_repo().read_pos.fetch_add(1, Ordering::Relaxed);
        // log::debug!("Add 1 to read_pos.");

        // if self.m_repo().read_pos.load(Ordering::Relaxed) as i32 >= PACK_NUM_LIMIT {
        //     self.m_repo().read_pos.store(0, Ordering::Relaxed);
        // }

        // lock.unlock();
        // mRepo.repoNotFull.notify_all();

        self.scan_pair_end(data)?;
        Ok(())
    }

    fn scan_pair_end(&self, pack: ReadPairPack) -> Result<bool, Box<dyn Error>> {
        let m_fusion_mapper = self.m_fusion_mapper_o.as_ref().unwrap();

        for (p, pair) in (0..(pack.count as usize)).zip(pack.data.into_iter()) {
            // let pair = pack.data.get(p).unwrap();
            let r1 = &pair.m_left;
            let r2 = &pair.m_right;

            let rcr1;
            let rcr2;

            let merged = pair.fast_merge();
            // if pair.m_left.m_name.contains(DBT) {
            //     log::debug!("merged={:#?}", merged);
            // };

            let mut mapable = false;

            let merged_rc;
            // if merged successfully, we only search the merged
            // log::debug!("p={}, merged={:?}", p, merged);
            if let Some(ref m) = merged {
                let mut match_merged = m_fusion_mapper.map_read(m, &mut mapable, 2, 20)?;
                // if pair.m_left.m_name.contains(DBT) {
                //     log::debug!("match_merged={:#?}", match_merged);
                // };

                // log::debug!("match_merged={:?}", match_merged);
                if let Some(mut mm) = match_merged {
                    mm.add_original_pair(pair.clone());
                    self.push_match(mm);
                } else if mapable {
                    merged_rc = m.reverse_complement();
                    let mut match_merged_rc =
                        m_fusion_mapper.map_read(&merged_rc, &mut mapable, 2, 20)?;
                    // if pair.m_left.m_name.contains(DBT) {
                    //     log::debug!("match_merged_rc={:#?}", match_merged_rc);
                    // };
                    if let Some(mut mmr) = match_merged_rc {
                        mmr.add_original_pair(pair.clone());
                        self.push_match(mmr);
                    }
                }
                continue;
            }
            // else still search R1 and R2 separatedly
            mapable = false;
            let match_r1 = m_fusion_mapper.map_read(r1, &mut mapable, 2, 20)?;
            // if pair.m_left.m_name.contains(DBT) {
            //     log::debug!("match_r1={:#?}", match_r1);
            // };
            if let Some(mut mr1) = match_r1 {
                mr1.add_original_pair(pair.clone());
                self.push_match(mr1);
            } else if mapable {
                rcr1 = r1.reverse_complement();
                let match_rcr1 = m_fusion_mapper.map_read(&rcr1, &mut mapable, 2, 20)?;
                // if pair.m_left.m_name.contains(DBT) {
                //     log::debug!("match_rcr1={:#?}", match_rcr1);
                // };
                if let Some(mut mrc1) = match_rcr1 {
                    mrc1.add_original_pair(pair.clone());
                    mrc1.set_reversed(true);
                    self.push_match(mrc1);
                }
            }

            mapable = false;

            let match_r2 = m_fusion_mapper.map_read(r2, &mut mapable, 2, 20)?;
            // if pair.m_left.m_name.contains(DBT) {
            //     log::debug!("match_r2={:#?}", match_r2);
            // };
            if let Some(mut mr2) = match_r2 {
                mr2.add_original_pair(pair.clone());
                self.push_match(mr2);
            } else if mapable {
                rcr2 = r2.reverse_complement();
                let match_rcr2 = m_fusion_mapper.map_read(&rcr2, &mut mapable, 2, 20)?;
                // if pair.m_left.m_name.contains(DBT) {
                //     log::debug!("match_rcr2={:#?}", match_rcr2);
                // };
                if let Some(mut mrc2) = match_rcr2 {
                    mrc2.add_original_pair(pair.clone());
                    mrc2.set_reversed(true);
                    self.push_match(mrc2);
                }
            }
        }

        Ok(true)
    }
    // #[track_caller]
    fn push_match(&self, m: ReadMatch) {
        // lock(self.m_fusion_mtx);

        // if m.m_read.m_name.contains(DBT) {
        // let loc = Location::caller();
        //     log::debug!(
        //         "push_match() called from {}:{}:{}",
        //         loc.file(),
        //         loc.line(),
        //         loc.column()
        //     );
        // }

        self.m_fusion_mapper_o.as_ref().unwrap().add_match(m);
        // lock.unlock();
    }

    pub(crate) fn text_report(&self) {
        todo!()
    }

    pub(crate) fn html_report(&mut self) -> Result<(), Box<dyn Error>> {
        if self.m_html_file == "" {
            return Ok(());
        }
        let mut reporter = HtmlReporter::new(
            self.m_html_file.clone(),
            self.m_fusion_mapper_o.as_mut().unwrap(),
        )?;

        reporter.run()?;

        Ok(())
    }

    pub(crate) fn json_report(&self) -> Result<(), Box<dyn Error>> {
        if self.m_json_file == "" {
            return Ok(());
        }
        let mut reporter = JsonReporter::new(
            self.m_json_file.clone(),
            self.m_fusion_mapper_o.as_ref().unwrap(),
        )?;

        reporter.run()?;

        Ok(())
    }
}

struct CPPMutex {}

#[cfg(test)]
mod test {
    #[test]
    fn cb_test() {
        use crossbeam::queue::ArrayQueue;

        let q = ArrayQueue::new(2);

        q.push('a');
        q.push('b');
        q.push('c');
        q.push('d');

        println!("{}", q.pop().unwrap());
    }
}
