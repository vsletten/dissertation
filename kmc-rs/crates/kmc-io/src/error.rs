//! The error type shared by every reader in this crate.
//!
//! The C++ model's entire error strategy is `Myerr::die(msg)` — print to
//! stderr and `exit(1)` from wherever the problem was noticed (`myerr.cpp`).
//! That conflates *detecting* a problem with *deciding what to do about it*.
//! The Rust convention separates the two: library code returns a value that
//! **is** the error, and only the binary decides whether that means exit,
//! retry, or a nicer message.

use std::fmt;
use std::path::PathBuf;

/// Anything that can go wrong while reading one of the `data.*` input files.
///
/// \[IDIOM\] Errors as enum values, not error codes or exceptions. A C++
/// function that "returns an int where -1 means failure" tells the caller
/// nothing and lets them forget to check. A `Result<T, ReadError>` return
/// type is a contract the compiler enforces: the caller *cannot* touch the
/// `T` without acknowledging the error case first (usually with `?`). The
/// enum also carries structured context — which file, what we were parsing —
/// so the message is assembled at the edge instead of formatted deep inside
/// the parser the way `Myerr::die` does.
#[derive(Debug)]
pub enum ReadError {
    /// The file could not be opened or read at the OS level.
    Io {
        /// Path we tried to read.
        path: PathBuf,
        /// The underlying OS error.
        ///
        /// \[IDIOM\] Error *chaining*: we wrap the lower-level error instead of
        /// flattening it to a string, so callers (and `{:#}`-style reporters)
        /// can walk the chain via [`std::error::Error::source`].
        source: std::io::Error,
    },
    /// The file opened fine but its contents did not parse.
    Parse {
        /// Path of the offending file.
        path: PathBuf,
        /// Human description of what was expected, mirroring the C++
        /// `Myerr::die` messages (e.g. "invalid number of steps in input
        /// file").
        what: String,
    },
}

impl fmt::Display for ReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ReadError::Io { path, source } => {
                // Mirrors futil.cpp: "could not open <name> for reading".
                write!(f, "could not read {}: {}", path.display(), source)
            }
            ReadError::Parse { path, what } => {
                write!(f, "{}: {}", path.display(), what)
            }
        }
    }
}

// [IDIOM] `std::error::Error` is a *trait* (an interface the type opts into),
// not a base class. C++ exception hierarchies force errors into an
// inheritance tree; Rust errors are plain values that additionally implement
// this trait so generic code (`Box<dyn Error>`, error reporters) can handle
// any of them uniformly. Note there is no `throw` anywhere — propagation is
// explicit dataflow through return values.
impl std::error::Error for ReadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            ReadError::Io { source, .. } => Some(source),
            ReadError::Parse { .. } => None,
        }
    }
}
