use std::{error::Error, fmt, io::Empty, panic::Location, str::FromStr, any::Any};

pub(crate) struct ErrorExplained {
    inner: Box<dyn Error>,
    loc: &'static Location<'static>,
}

impl Error for ErrorExplained {}

impl fmt::Display for ErrorExplained {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self.inner, f).and_then(|_| {
            let loc = self.loc;
            write!(f, " at {}:{}:{}", loc.file(), loc.line(), loc.column())
        })
    }
}

impl fmt::Debug for ErrorExplained {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Debug::fmt(&self.inner, f).and_then(|_| {
            let loc = self.loc;
            write!(f, " at {}:{}:{}", loc.file(), loc.line(), loc.column())
        })
    }
}

pub(crate) trait OrExaplain<T> {
    /// ok or explain err
    fn or_exp(self) -> Result<T, ErrorExplained>;
}

impl<T, E> OrExaplain<T> for Result<T, E>
where
    E: Error + 'static,
{
    #[track_caller]
    fn or_exp(self) -> Result<T, ErrorExplained> {
        match self {
            Ok(v) => Ok(v),
            Err(err) => Err(ErrorExplained {
                inner: From::from(err),
                loc: Location::caller(),
            }),
        }
    }
}




// impl<E:Error> From<E> for ErrorExplained
// {
//     #[track_caller]
//     fn from(value: E) -> Self {
//         ErrorExplained {
//             inner: value,
//             loc: Location::caller(),
//         }
//     }
// }

macro_rules! make_custom_error {
    ($n:ident, $error_msg:tt) => {
        pub(crate) struct $n {}

        impl $n {
            pub(crate) fn new() -> Self {
                Self {}
            }
        }

        impl fmt::Display for $n {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, $error_msg)
            }
        }

        impl fmt::Debug for $n {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, "{}: {}", stringify!($n), $error_msg)
            }
        }

        impl Error for $n {}
    };
}



pub(crate) use make_custom_error;
macro_rules! make_custom_error2 {
    ($n:ident, $($error_msg:tt)*) => {
        pub(crate) struct $n<C> {
            cxt: C,
        }

        impl<C:fmt::Debug> $n<C> {
            pub(crate) fn new(cxt: C) -> Self {
                Self {
                    cxt
                }
            }
        }

        impl<C:fmt::Debug> fmt::Display for $n<C> {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, $($error_msg)*)
                    // .and_then(|_| write!(f, " cxt: "))
                    // .and_then(|_| fmt::Display::fmt(&self.cxt, f))
            }
        }

        impl<C:fmt::Debug> fmt::Debug for $n<C> {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.debug_struct(stringify!($n)).field("cxt", &self.cxt).finish()
            }
            // fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            //     write!(f, "{}: ", stringify!($n))
            //         .and_then(|_| fmt::Display::fmt(self, f))
                    
            // }
        }

        impl<C:fmt::Debug> Error for $n<C> {}
    };
}

pub(crate) use make_custom_error2;

// pub(crate) trait DebugExceptString
// where
//     Self: fmt::Debug,
// {
//     fn des(&self) -> String {
//         format!("{:?}", self)
//     }
// }

// impl DebugExceptString for String {
//     fn des(&self) -> String {
//         format!("{}", self)
//     }
// }

// impl<T> DebugExceptString for T {}

// fn debug_except_string<T: Any + fmt::Debug>(v: &T) {
//     let v= v as &dyn Any;
//     match v.downcast_ref::<String>() {
//         Some(a) => a.clone(),
//         None => {
//             match v.downcast_ref::<str>() {
//                 Some(a) => a.to_owned(),
//                 None => todo!(),
//             }
//         },
//     }
// }

macro_rules! make_custom_error3 {
    ($n:ident, $($error_msg:tt)*) => {
        pub(crate) struct $n {
            cxt: String,
        }

        impl $n {
            pub(crate) fn new<C: fmt::Debug>(cxt: &C) -> Self {
                Self {
                    cxt: format!("{:?}", cxt),
                }
            }
        }

        impl fmt::Display for $n {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, $($error_msg)*)
                    // .and_then(|_| write!(f, " cxt: "))
                    // .and_then(|_| fmt::Display::fmt(&self.cxt, f))
            }
        }

        impl fmt::Debug for $n {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.debug_struct(stringify!($n)).field("cxt", &self.cxt).finish()
            }
            // fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            //     write!(f, "{}: ", stringify!($n))
            //         .and_then(|_| fmt::Display::fmt(self, f))
                    
            // }
        }

        impl Error for $n {}
    };
}

pub(crate) use make_custom_error3;


macro_rules! make_custom_error4 {
    ($n:ident, $($error_msg:tt)*) => {
        pub(crate) struct $n {
            cxt: String,
            loc: String,
        }

        impl $n {
            #[track_caller]
            pub(crate) fn new<C: fmt::Debug>(cxt: &C) -> Self {
                let loc = Location::caller();
                Self {
                    cxt: format!("{:?}", cxt),
                    loc: format!("{}:{}:{}", loc.file(), loc.line(), loc.column()),
                }
            }
        }

        impl fmt::Display for $n {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, $($error_msg)*)
                    // .and_then(|_| write!(f, " {}:{}:{}", ))
                    // .and_then(|_| fmt::Display::fmt(&self.cxt, f))
            }
        }

        impl fmt::Debug for $n {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.debug_struct(stringify!($n)).field("cxt", &self.cxt).field("loc", &self.loc).finish()
            }
            // fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            //     write!(f, "{}: ", stringify!($n))
            //         .and_then(|_| fmt::Display::fmt(self, f))
                    
            // }
        }

        impl Error for $n {}
    };
}

pub(crate) use make_custom_error4;




#[cfg(test)]
mod test {
    use std::{error::Error, io};
    use std::fmt;

    use super::*;

    fn t2() -> Result<i32, ErrorExplained> {
        Ok("1.023".parse::<i32>().or_exp()?)
    }
    #[test]
    fn frome() {
        t2().unwrap();
    }

    fn _fromerr() -> Result<(), Box<dyn Error>> {
        let mut s = "AAA".to_string();

        make_custom_error3!(Error1, "this is error1.");

        Err(Error1::new(&mut s))?


    }

    #[test]
    fn test_ce_macro() {
        _fromerr().unwrap();

        make_custom_error2!(Error1, "this is error1.");

        let e = Error1::new(1);

        let e2 = io::Error::new(io::ErrorKind::AlreadyExists, "Oh HO!");

        println!("{e2}", );
        println!("{e2:?}", );
        println!("{e}");
        println!("{e:?}");
    }
}
