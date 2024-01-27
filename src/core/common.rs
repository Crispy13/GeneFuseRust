use super::gene::Gene;

#[derive(Clone, PartialEq, PartialOrd, Debug)]
pub(crate) struct GenePos {
    pub(crate) contig: i16,
    pub(crate) position: i32,
}

impl Default for GenePos {
    fn default() -> Self {
        Self {
            contig: 0,
            position: 0,
        }
    }
}

// the limit of the queue to store the packs
// error may happen if it generates more packs than this number
pub(crate) const PACK_NUM_LIMIT: i32 = 5000000;

// how many reads one pack has
pub(crate) const PACK_SIZE: i32 = 1000;

// if one pack is produced, but not consumed, it will be kept in the memory
// this number limit the number of in memory packs
// if the number of in memory packs is full, the producer thread should sleep
pub(crate) const PACK_IN_MEM_LIMIT: i32 = 100;

// the key dup in normal level will be kept, in high level will be skipped
pub(crate) const DUPE_NORMAL_LEVEL: i16 = -1;
pub(crate) const DUPE_HIGH_LEVEL: i16 = -2;
