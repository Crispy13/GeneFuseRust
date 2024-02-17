use std::{
    error,
    sync::{atomic::{AtomicBool, AtomicUsize, Ordering}, Arc},
    thread::sleep,
    time::Duration,
};

use crossbeam::queue::ArrayQueue;
use rayon::{ThreadPool, ThreadPoolBuilder};

use crate::core::{
    common::{PACK_NUM_LIMIT, PACK_SIZE},
    fastq_reader::FastqReader,
    sequence::reverse_complement,
};

use super::{
    fasta_reader::FastaReader, fusion_mapper::FusionMapper, fusion_scan::Error, html_reporter::HtmlReporter, json_reporter::JsonReporter, read::SequenceRead, read_match::ReadMatch
};

#[derive(Debug)]
struct ReadPack {
    data: Vec<SequenceRead>,
    count: i32,
}

struct ReadRepository {
    pack_buffer: ArrayQueue<ReadPack>,
    read_pos: AtomicUsize,
    write_pos: AtomicUsize,
    read_counter: AtomicUsize,
    // mtx: CPPMutex,
    // read_counter_mtx: CPPMutex,
    // repo_not_full: Condvar,
    // repo_not_empty: Condvar,
}

pub(crate) struct SingleEndScanner {
    m_fusion_file: String,
    m_ref_file: String,
    m_read1_file: String,
    m_html_file: String,
    m_json_file: String,
    m_repo_o: Option<ReadRepository>,
    m_produce_finished: AtomicBool,
    m_thread_num: i32,
    m_fusion_mapper_o: Option<FusionMapper>,
    m_thread_pool: ThreadPool,
}

impl SingleEndScanner {
    pub(crate) fn new(
        fusion_file: String,
        ref_file: String,
        read1_file: String,
        html: String,
        json: String,
        thread_num: i32,
    ) -> Self {
        let itp = ThreadPoolBuilder::new().num_threads(thread_num as usize).build().unwrap();

        Self {
            m_fusion_file: fusion_file,
            m_ref_file: ref_file,
            m_read1_file: read1_file,
            m_html_file: html,
            m_json_file: json,
            m_repo_o: None,
            m_produce_finished: AtomicBool::new(false),
            m_thread_num: thread_num,
            m_fusion_mapper_o: None,
            m_thread_pool: itp,
            // repo_not_full: Condvar::new(),
            // repo_not_empty: Condvar::new(),
        }
    }

    fn m_repo(&self) -> &ReadRepository {
        self.m_repo_o.as_ref().unwrap()
    }

    fn _scan(&mut self) -> Result<bool, Error> {
        self.m_thread_pool.scope(|tps| {
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
        m_fusion_mapper.filter_matches(&self.m_thread_pool);
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

    pub(crate) fn scan(&mut self) -> Result<bool, Error> {
        log::debug!("Entered into scan.");
        self.m_fusion_mapper_o = Some(FusionMapper::from_ref_and_fusion_files(
            &self.m_ref_file,
            &self.m_fusion_file,
        )?);
        log::debug!("Made fusion mapper.");

        self.init_pack_repository();

        Ok(self._scan()?)
    }

    fn push_match(&self, m: ReadMatch) {
        self.m_fusion_mapper_o.as_ref().unwrap().add_match(m);
    }

    fn scan_single_end(&self, pack: ReadPack) -> Result<bool, Error> {
        let m_fusion_mapper = self.m_fusion_mapper_o.as_ref().unwrap();

        for (p, r1) in (0..(pack.count as usize)).zip(pack.data.into_iter()) {
            let mut mapable = false;
            let match_r1 = m_fusion_mapper.map_read(&r1, &mut mapable, 2, 20)?;

            if let Some(mut mr1) = match_r1 {
                mr1.add_original_read(r1.clone());
                self.push_match(mr1);
            } else if mapable {
                let rcr1 = r1.reverse_complement();
                let match_rcr1 = m_fusion_mapper.map_read(&rcr1, &mut mapable, 2, 20)?;
                if let Some(mut mrcr1) = match_rcr1 {
                    mrcr1.add_original_read(r1.clone());
                    mrcr1.set_reversed(true);
                    self.push_match(mrcr1);
                }
            }
        }

        Ok(true)
    }

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
            let repo = ReadRepository {
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

    fn produce_pack(&self, mut pack: ReadPack) -> Result<(), Error> {
        log::debug!("Entered produce_pack()");

        let m_repo = self.m_repo_o.as_ref().unwrap();

        loop {
            match m_repo.pack_buffer.push(pack) {
                Ok(_) => break Ok(()),
                Err(p) => {
                    if self.m_thread_num == 1 {
                        self.consume_pack()?;
                    }
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

    fn consume_pack(&self) -> Result<(), Error> {
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

        self.scan_single_end(data)?;
        Ok(())
    }

    fn producer_task(&self) -> Result<(), Error> {
        log::debug!("Entered producer_task()");

        let mut slept = 0;
        let mut data = Vec::<SequenceRead>::with_capacity(PACK_SIZE as usize);

        let mut reader = FastqReader::new(&self.m_read1_file, true)?;

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
                let pack = ReadPack { data, count };
                self.produce_pack(pack)?;

                //re-initialize data for next pack
                data = Vec::<SequenceRead>::with_capacity(PACK_SIZE as usize);

                // reset count to 0
                count = 0;
                // if the consumer is far behind this producer, sleep and wait to limit memory usage
                // while self.m_repo().write_pos.load(Ordering::Relaxed) as i32
                //     - (self.m_repo().read_pos.load(Ordering::Relaxed) as i32)
                //     > PACK_NUM_LIMIT
                // {
                //     slept += 1;
                //     sleep(Duration::from_micros(1000));
                // }
            }
        }

        let pack = ReadPack { data, count };
        self.produce_pack(pack)?;

        // lock self.m_repo().read_counter_mtx;
        log::debug!("producing finished.");
        self.m_produce_finished.store(true, Ordering::Relaxed);
        // lock.unlock();

        // producing tasks has been done, from now try to consume the tasks.
        while !self.m_repo_o.as_ref().unwrap().pack_buffer.is_empty() {
            self.consume_pack()?;
        }

        Ok(())
    }

    fn consumer_task(&self) -> Result<(), Error> {
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

    pub(crate) fn text_report(&self) {
        todo!()
    }

    pub(crate) fn html_report(&mut self) -> Result<(), Error> {
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

    pub(crate) fn json_report(&self) -> Result<(), Error> {
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

    pub(crate) fn scan_per_fusion_csv(
        &mut self,
        fasta_reader: Arc<FastaReader>,
    ) -> Result<bool, Error> {
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

    // pub(crate) fn drop_and_get_back_fasta_reader(self) -> FastaReader {
    //     self.m_fusion_mapper_o
    //         .unwrap()
    //         .m_indexer
    //         .m_reference
    //         .unwrap()
    // }
}
