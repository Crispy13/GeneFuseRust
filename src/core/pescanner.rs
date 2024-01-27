use rayon::{prelude::*, ThreadPoolBuilder};
use std::{
    error::Error,
    sync::{
        atomic::{AtomicBool, AtomicUsize, Ordering},
        Condvar, Mutex, RwLock,
    },
    thread::sleep,
    time::Duration,
};

use super::{
    common::{PACK_NUM_LIMIT, PACK_SIZE},
    fastq_reader::FastqReaderPair,
    fusion_mapper::FusionMapper,
    read::SequenceReadPair,
    read_match::ReadMatch,
};
use crate::core::{html_reporter::HtmlReporter, json_reporter::JsonReporter};
use crossbeam::queue::ArrayQueue;

#[derive(Debug)]
struct ReadPairPack {
    data: Vec<SequenceReadPair>,
    count: i32,
}

struct CPPConditionVariable {}

struct ReadPairRepository {
    pack_buffer: ArrayQueue<ReadPairPack>,
    read_pos: AtomicUsize,
    write_pos: AtomicUsize,
    read_counter: AtomicUsize,
    mtx: CPPMutex,
    read_counter_mtx: CPPMutex,
    repo_not_full: Condvar,
    repo_not_empty: Condvar,
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
    m_fusion_mtx: CPPMutex,
    m_thread_num: i32,
    m_fusion_mapper_o: Option<Mutex<FusionMapper>>,
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
            m_fusion_mtx: CPPMutex {},
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
                mtx: CPPMutex {},
                read_counter_mtx: CPPMutex {},
                repo_not_full: Condvar::new(),
                repo_not_empty: Condvar::new(),
            };

            self.m_repo_o = Some(repo);
        }
    }

    fn produce_pack(&self, pack: ReadPairPack) {
        log::debug!("Entered produce_pack()");
        // lock m_repo.mtx;

        // while (self.m_repo().write_pos + 1) % PACK_NUM_LIMIT as usize == self.m_repo().read_pos {
        //     // self.m_repo().repo_not_full.wait(lock);
        // }

        // let wp = self.m_repo().write_pos;
        let m_repo = self.m_repo_o.as_ref().unwrap();

        while m_repo.pack_buffer.is_full() {
            // m_repo.repo_not_full.wait(m_repo);
        }

        m_repo.pack_buffer.push(pack).unwrap(); // *self.m_repo_mut().pack_buffer.get_mut(wp).unwrap() = pack;
        m_repo.write_pos.fetch_add(1, Ordering::Relaxed);
        log::debug!("Add 1 to write_pos.");

        if m_repo.write_pos.load(Ordering::Relaxed) == PACK_NUM_LIMIT as usize {
            m_repo.write_pos.store(0, Ordering::Relaxed);
        }

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
                // if the consumer is far behind this producer, sleep and wait to limit memory usage
                while self.m_repo().write_pos.load(Ordering::Relaxed) as i32
                    - (self.m_repo().read_pos.load(Ordering::Relaxed) as i32)
                    > PACK_NUM_LIMIT
                {
                    slept += 1;
                    sleep(Duration::from_micros(1000));
                }
            }
        }

        let pack = ReadPairPack { data, count };
        self.produce_pack(pack);

        // lock self.m_repo().read_counter_mtx;
        self.m_produce_finished.store(true, Ordering::Relaxed);
        // lock.unlock();

        Ok(())
    }

    pub(crate) fn scan(&mut self) -> Result<bool, Box<dyn Error>> {
        log::debug!("Entered into scan.");
        self.m_fusion_mapper_o = Some(Mutex::new(FusionMapper::from_ref_and_fusion_files(
            &self.m_ref_file,
            &self.m_fusion_file,
        )?));
        log::debug!("Made fusion mapper.");

        self.init_pack_repository();

        let tp = ThreadPoolBuilder::new()
            .num_threads(self.m_thread_num as usize)
            .build()
            .unwrap();

        tp.scope(|tps| {
            tps.spawn(|tps| self.producer_task().unwrap());

            for t in (0..(tp.current_num_threads() - 1)) {
                tps.spawn(|tps| {
                    self.consumer_task().unwrap();
                })
            }
        });

        let mut m_fusion_mapper = self.m_fusion_mapper_o.as_ref().unwrap().lock().unwrap();

        m_fusion_mapper.filter_matches();
        m_fusion_mapper.sort_matches();
        m_fusion_mapper.cluster_matches();

        self.html_report().unwrap();
        self.json_report().unwrap();

        self.m_fusion_mapper_o
            .as_ref()
            .unwrap()
            .lock()
            .unwrap()
            .free_matches();

        Ok(true)
    }

    fn consumer_task(&self) -> Result<(), Box<dyn Error>> {
        log::debug!("Entered consumer_task()");
        loop {
            // lock = self.m_repo().read_counter_mtx;
            if self.m_produce_finished.load(Ordering::Relaxed)
                && self.m_repo().write_pos.load(Ordering::Relaxed)
                    == self.m_repo().read_pos.load(Ordering::Relaxed)
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

        Ok(())
    }

    fn consume_pack(&self) -> Result<(), Box<dyn Error>> {
        // let data = ReadPairPack {
        //     data: todo!(),
        //     count: todo!(),
        // };

        // lock self.m_repo().mtx;

        // read buffer is empty, just wait here.
        while self.m_repo().write_pos.load(Ordering::Relaxed) as i32 % PACK_NUM_LIMIT
            == self.m_repo().read_pos.load(Ordering::Relaxed) as i32 % PACK_NUM_LIMIT
        {
            if self.m_produce_finished.load(Ordering::Relaxed) {
                // lock.unlock();
                return Ok(());
            }
            // self.m_repo().repo_not_empty.wait(lock);
        }

        // let data = self
        //     .m_repo_o
        //     .as_ref()
        //     .unwrap()
        //     .pack_buffer
        //     .get(self.m_repo().read_pos)
        //     .unwrap();

        while self.m_repo_o.as_ref().unwrap().pack_buffer.is_empty() {}

        let data = self.m_repo_o.as_ref().unwrap().pack_buffer.pop().unwrap();

        self.m_repo().read_pos.fetch_add(1, Ordering::Relaxed);
        log::debug!("Add 1 to read_pos.");

        if self.m_repo().read_pos.load(Ordering::Relaxed) as i32 >= PACK_NUM_LIMIT {
            self.m_repo().read_pos.store(0, Ordering::Relaxed);
        }

        // lock.unlock();
        // mRepo.repoNotFull.notify_all();

        self.scan_pair_end(data)?;
        Ok(())
    }

    fn scan_pair_end(&self, pack: ReadPairPack) -> Result<bool, Box<dyn Error>> {
        let mut m_fusion_mapper = self.m_fusion_mapper_o.as_ref().unwrap().lock().unwrap();

        for (p, pair) in (0..(pack.count as usize)).zip(pack.data.into_iter()) {
            // let pair = pack.data.get(p).unwrap();
            let r1 = &pair.m_left;
            let r2 = &pair.m_right;

            let rcr1;
            let rcr2;

            let merged = pair.fast_merge();
            let mut mapable = false;

            let merged_rc;
            // if merged successfully, we only search the merged
            if let Some(ref m) = merged {
                let mut match_merged = m_fusion_mapper.map_read(m, &mut mapable, 2, 20)?;
                if let Some(mut mm) = match_merged {
                    mm.add_original_pair(pair.clone());
                    self.push_match(mm);
                } else if mapable {
                    merged_rc = m.reverse_complement();
                    let mut match_merged_rc =
                        m_fusion_mapper.map_read(&merged_rc, &mut mapable, 2, 20)?;

                    if let Some(mut mmr) = match_merged_rc {
                        mmr.add_original_pair(pair.clone());
                        self.push_match(mmr);
                    }
                }
                continue;
            }
            // else still search R1 and R2 separatedly
            mapable = false;
            let mut match_r1 = m_fusion_mapper.map_read(r1, &mut mapable, 2, 20)?;
            if let Some(mut mr1) = match_r1 {
                mr1.add_original_pair(pair.clone());
                self.push_match(mr1);
            } else if mapable {
                rcr1 = r1.reverse_complement();
                let match_rcr1 = m_fusion_mapper.map_read(&rcr1, &mut mapable, 2, 20)?;
                if let Some(mut mrc1) = match_rcr1 {
                    mrc1.add_original_pair(pair.clone());
                    mrc1.set_reversed(true);
                    self.push_match(mrc1);
                }
            }

            mapable = false;

            let match_r2 = m_fusion_mapper.map_read(r2, &mut mapable, 2, 20)?;
            if let Some(mut mr2) = match_r2 {
                mr2.add_original_pair(pair.clone());
                self.push_match(mr2);
            } else if mapable {
                rcr2 = r2.reverse_complement();
                let match_rcr2 = m_fusion_mapper.map_read(&rcr2, &mut mapable, 2, 20)?;
                if let Some(mut mrc2) = match_rcr2 {
                    mrc2.add_original_pair(pair.clone());
                    mrc2.set_reversed(true);
                    self.push_match(mrc2);
                }
            }
        }

        Ok(true)
    }
    fn push_match(&self, m: ReadMatch) {
        // lock(self.m_fusion_mtx);
        self.m_fusion_mapper_o
            .as_ref()
            .unwrap()
            .lock()
            .unwrap()
            .add_match(m.clone());
        // lock.unlock();
    }

    pub(crate) fn text_report(&self) {
        todo!()
    }

    pub(crate) fn html_report(&self) -> Result<(), Box<dyn Error>> {
        if self.m_html_file == "" {
            return Ok(());
        }
        let mut reporter = HtmlReporter::new(
            self.m_html_file.clone(),
            self.m_fusion_mapper_o.as_ref().unwrap().lock().unwrap(),
        )?;

        reporter.run();

        Ok(())
    }

    pub(crate) fn json_report(&self) -> Result<(), Box<dyn Error>> {
        if self.m_json_file == "" {
            return Ok(());
        }
        let mut reporter = JsonReporter::new(
            self.m_json_file.clone(),
            self.m_fusion_mapper_o.as_ref().unwrap().lock().unwrap(),
        )?;

        reporter.run();

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
