use std::io::{Take, BufReader, Read, BufRead};
use std::fmt;
use std::error::Error;


use super::he::{make_custom_error3};

make_custom_error3!(InsufficientTakeAmountError, "Insufficient take amount.");

pub(crate) struct LimitedBufReader<T> {
    inner: Take<T>,
    limit: u64,
}

impl<T:BufRead> LimitedBufReader<T> {
    pub(crate) fn new(inner: T, limit: u64) -> Self {
        Self {
            inner: inner.take(limit),
            limit: limit,
        }
    }
}

impl<T:BufRead> Read for LimitedBufReader<T> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        self.inner.read(buf)
    }
}

impl<T:BufRead> BufRead for LimitedBufReader<T> {
    fn read_until(&mut self, byte: u8, buf: &mut Vec<u8>) -> std::io::Result<usize> {
        match self.inner.read_until(byte, buf) {
            Ok(rl) => {
                // Check unexected truncating happened.
                if *buf.last().unwrap() != b'\n' {
                    self.inner.set_limit(1);
                    match self.inner.fill_buf() {
                        Ok(bs) => {if bs.len() > 0 {
                            panic!("{:?}", InsufficientTakeAmountError::new(&(buf, *bs.first().unwrap() as char)));
                        }},
                        Err(err) => Err(err)?,
                    }
                }

                self.inner.set_limit(self.limit); // reset limit to read more than once.
                Ok(rl)
            },
            Err(err) => {
                Err(err)
            },
        }
    }

    fn read_line(&mut self, buf: &mut String) -> std::io::Result<usize> {
        match self.inner.read_line(buf) {
            Ok(rl) => {
                // Check unexected truncating happened.
                if *buf.as_bytes().last().unwrap() != b'\n' {
                    self.inner.set_limit(1);
                    match self.inner.fill_buf() {
                        Ok(bs) => {if bs.len() > 0 {
                            panic!("{:?}", InsufficientTakeAmountError::new(&(buf, *bs.first().unwrap() as char)));
                        }},
                        Err(err) => Err(err)?,
                    }
                }

                self.inner.set_limit(self.limit); // reset limit to read more than once.
                Ok(rl)
            },
            Err(err) => {
                Err(err)
            },
        }
    }

    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        self.inner.fill_buf()
    }

    fn consume(&mut self, amt: usize) {
        self.inner.consume(amt)
    }
}
