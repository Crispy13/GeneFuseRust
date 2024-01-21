pub(crate) struct Array2D<T, const N:usize> {
    inner: [T; N],
    shape: (usize, usize),
}

impl<T, const N:usize> Array2D<T, N> {
    pub(crate) fn get(&self, r:usize, c:usize) -> Option<&T> {
        self.inner.get(Self::get_flatten_idx(&self.shape, r, c))
    }

    pub(crate) fn get_mut(&mut self, r:usize, c:usize) -> Option<&mut T> {
        self.inner.get_mut(Self::get_flatten_idx(&self.shape, r, c))
    }

    #[inline]
    fn get_flatten_idx(shape:&(usize, usize), r:usize, c:usize) -> usize {
        r * shape.1 + c
    }
}

