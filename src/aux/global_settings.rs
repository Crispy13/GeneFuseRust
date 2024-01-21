use std::sync::{OnceLock, RwLock, RwLockReadGuard, RwLockWriteGuard};

pub(crate) struct GlobalSettings {
    pub(crate) marked_only_for_vcf: bool,
    pub(crate) unique_requirement: usize,
    pub(crate) deletion_threshold: usize,
    pub(crate) output_deletions: bool,
    pub(crate) output_untranslated: bool,
    pub(crate) skip_key_dup_threshold: usize,
    pub(crate) major_gene_key_requirement: i32,
    pub(crate) minor_gene_key_requirement: i32,
    pub(crate) mismatch_threshold: i32,
}

impl Default for GlobalSettings {
    fn default() -> Self {
        Self {
            marked_only_for_vcf: false,
            unique_requirement: 2,
            deletion_threshold: 50,
            output_deletions: false,
            output_untranslated: false,
            skip_key_dup_threshold: 5,
            major_gene_key_requirement: 40,
            minor_gene_key_requirement: 20,
            mismatch_threshold: 10,
        }
    }
}

impl GlobalSettings {
    #[inline]
    pub(crate) fn set_marked_only_for_vcf(&mut self, flag: bool) {
        self.marked_only_for_vcf = flag;
    }

    #[inline]
    pub(crate) fn set_unique_requirement(&mut self, val: usize) {
        self.unique_requirement = val;
    }

    #[inline]
    pub(crate) fn set_deletion_threshold(&mut self, val: usize) {
        self.deletion_threshold = val;
    }

    #[inline]
    pub(crate) fn set_output_deletions(&mut self, flag: bool) {
        self.output_deletions = flag;
    }

    #[inline]
    pub(crate) fn set_output_untranslated(&mut self, flag: bool) {
        self.output_untranslated = flag;
    }
}

static GLOBAL_SETTINGS: OnceLock<RwLock<GlobalSettings>> = OnceLock::new();
pub(crate) fn global_settings() -> RwLockReadGuard<'static, GlobalSettings> {
    GLOBAL_SETTINGS.get_or_init(|| RwLock::new(GlobalSettings::default())).read().unwrap()
}

pub(crate) fn global_settings_w() -> RwLockWriteGuard<'static, GlobalSettings> {
    GLOBAL_SETTINGS.get_or_init(|| RwLock::new(GlobalSettings::default())).write().unwrap()
}
