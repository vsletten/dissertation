//! Tokenizer for the legacy input-file dialect.
//!
//! The C++ readers are built from exactly two primitives, and this module
//! ports both **with the same division of labor** so every reader can be a
//! line-by-line transliteration of its C++ counterpart:
//!
//! 1. `std::ifstream >> value` — skip whitespace, read one token, convert.
//!    Ported as [`Scanner::next_i32`] / [`Scanner::next_f32`].
//! 2. `Futil::EatComment(f, '#')` — skip whitespace; if the next character
//!    is `#`, discard to end of line; repeat (handles stacked comment
//!    lines). Ported as [`Scanner::eat_comment`].
//!
//! **Faithfulness note:** it would be *simpler* to strip `#` comments
//! globally and split on whitespace. We deliberately do not. The C++ only
//! tolerates comments at the exact points where `EatComment` is called — a
//! comment anywhere else makes `operator>>` fail and the program die. A
//! global-stripping scanner would accept files the C++ rejects, which is the
//! kind of quiet behavioral drift the faithful spec exists to prevent. So
//! `next_*` refuses to skip comments, and each reader calls `eat_comment` at
//! precisely the C++ call sites.

use std::path::{Path, PathBuf};

use crate::error::ReadError;

/// A cursor over the full text of one input file.
///
/// The C++ streams the file; we read it whole (these files are a few KB) and
/// walk a byte index. Same observable behavior, much simpler lifetime story.
pub struct Scanner {
    /// Raw file bytes. The legacy files are ASCII; scanning bytes instead of
    /// `char`s sidesteps UTF-8 boundary questions entirely.
    text: Vec<u8>,
    /// Current position in `text`.
    pos: usize,
    /// Which file we're reading — carried so every error names its source.
    path: PathBuf,
}

impl Scanner {
    /// Open `path` and return a scanner over its contents.
    ///
    /// Dies with `ReadError::Io` if the file can't be read — the moral
    /// equivalent of `Futil::OpenInputFile`'s `Myerr::die`, except the
    /// *caller* chooses whether that is fatal.
    pub fn open(path: &Path) -> Result<Self, ReadError> {
        // [IDIOM] `map_err` + `?`: convert the library-level error into our
        // domain error, then propagate. This two-step is the Rust spelling of
        // "catch, wrap, rethrow" — but it's visible in the signature and
        // costs no hidden control flow: `?` is an early `return Err(...)`,
        // nothing more. When reviewing agent-written Rust, grep for
        // `.unwrap()` where a `?` belongs — unwrap converts a recoverable
        // error into a crash and is the #1 shortcut agents take.
        let text = std::fs::read(path).map_err(|source| ReadError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        Ok(Scanner {
            text,
            pos: 0,
            path: path.to_path_buf(),
        })
    }

    /// Build a scanner from in-memory text (for tests).
    ///
    /// `label` stands in for the file name in error messages.
    ///
    /// \[IDIOM\] `impl Into<PathBuf>` — generic over "anything convertible to a
    /// PathBuf" (`&str`, `String`, `PathBuf`...). This is compile-time
    /// polymorphism like a C++ template parameter, but constrained by a
    /// declared trait bound instead of duck typing: the function body can
    /// only use what `Into<PathBuf>` promises, so misuse fails at the call
    /// site with a readable error, not three template layers deep.
    pub fn from_text(text: &str, label: impl Into<PathBuf>) -> Self {
        Scanner {
            text: text.as_bytes().to_vec(),
            pos: 0,
            path: label.into(),
        }
    }

    /// True if `b` is whitespace by C's `isspace` (space, \t, \n, \v, \f, \r)
    /// — the exact set `operator>>` and `EatComment` skip.
    fn is_space(b: u8) -> bool {
        matches!(b, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
    }

    /// Advance past any whitespace.
    fn skip_ws(&mut self) {
        while self.pos < self.text.len() && Self::is_space(self.text[self.pos]) {
            self.pos += 1;
        }
    }

    /// Port of `Futil::EatComment(f, '#')`: skip whitespace, and if the next
    /// character starts a `#` comment, discard through end of line; loop so
    /// stacked comment lines (and blank lines between them) all disappear.
    /// If the next non-space character is *not* `#`, it is left in place
    /// (the C++ `putback`).
    pub fn eat_comment(&mut self) {
        loop {
            self.skip_ws();
            if self.pos < self.text.len() && self.text[self.pos] == b'#' {
                while self.pos < self.text.len() && self.text[self.pos] != b'\n' {
                    self.pos += 1;
                }
                // loop again: there may be another comment line
            } else {
                return;
            }
        }
    }

    /// Read the next whitespace-delimited token, or `None` at end of input.
    ///
    /// \[IDIOM\] `Option<&str>` for "might not be there". C++ signals stream
    /// exhaustion through the stream's fail state — an out-of-band flag you
    /// must remember to test (`if (!(f >> x)) die(...)`). `Option` moves that
    /// flag *into the type*: there is no way to use the token without first
    /// deciding what `None` means here.
    fn next_token(&mut self) -> Option<&str> {
        self.skip_ws();
        let start = self.pos;
        while self.pos < self.text.len() && !Self::is_space(self.text[self.pos]) {
            self.pos += 1;
        }
        if self.pos == start {
            None
        } else {
            // The input dialect is ASCII; from_utf8 on a byte slice of it
            // cannot fail, but we still avoid unsafe/unchecked conversion.
            std::str::from_utf8(&self.text[start..self.pos]).ok()
        }
    }

    /// Read one token and parse it as `i32` — the port of `f >> some_int`.
    ///
    /// `what` mirrors the C++ die-messages ("invalid number of steps in
    /// input file") so diagnostics stay recognizable to anyone who knew the
    /// old model.
    pub fn next_i32(&mut self, what: &str) -> Result<i32, ReadError> {
        let path = self.path.clone();
        let tok = self.next_token();
        match tok.and_then(|t| t.parse::<i32>().ok()) {
            Some(v) => Ok(v),
            None => Err(ReadError::Parse {
                path,
                what: what.to_string(),
            }),
        }
    }

    /// Read one token and parse it as `f32` — the port of `f >> some_float`.
    ///
    /// Numeric fidelity: Rust's `str::parse::<f32>` and libstdc++'s
    /// `num_get<float>` both produce the **correctly rounded** nearest f32
    /// for a decimal token, so parsed constants are bit-identical between
    /// the two implementations. (This matters: the M3 bitwise gate rides on
    /// every parsed coordinate.)
    pub fn next_f32(&mut self, what: &str) -> Result<f32, ReadError> {
        let path = self.path.clone();
        let tok = self.next_token();
        match tok.and_then(|t| t.parse::<f32>().ok()) {
            Some(v) => Ok(v),
            None => Err(ReadError::Parse {
                path,
                what: what.to_string(),
            }),
        }
    }
}

// [IDIOM] `#[cfg(test)] mod tests` — unit tests live in the same file as the
// code, compiled only for `cargo test`. Coming from C++ (separate test
// binaries, external frameworks), the important property is cultural: tests
// are zero-setup, so "every parser change lands with a test" has no excuse
// not to happen. Reviewer's note: a Rust module with no `#[cfg(test)]` block
// at the bottom is the first thing to ask an agent about.
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tokens_and_comments_interleave_like_the_cpp() {
        let mut s = Scanner::from_text("42  # trailing comment\n# full line\n  7 rest", "test");
        assert_eq!(s.next_i32("first").unwrap(), 42);
        s.eat_comment(); // eats BOTH comment lines, like EatComment's loop
        assert_eq!(s.next_i32("second").unwrap(), 7);
    }

    #[test]
    fn comments_are_not_skipped_without_eat_comment() {
        // Faithfulness: without EatComment, the C++ would try to parse '#'
        // as a number and die. Our next_i32 must fail the same way.
        let mut s = Scanner::from_text("1 # comment\n2", "test");
        assert_eq!(s.next_i32("one").unwrap(), 1);
        assert!(s.next_i32("should hit the comment").is_err());
    }

    #[test]
    fn eof_yields_parse_error_not_panic() {
        let mut s = Scanner::from_text("  \n\t", "test");
        assert!(s.next_i32("nothing left").is_err());
    }

    #[test]
    fn floats_parse_correctly_rounded() {
        let mut s = Scanner::from_text("-0.0262374 8.92893", "test");
        assert_eq!(s.next_f32("a").unwrap(), -0.0262374_f32);
        assert_eq!(s.next_f32("b").unwrap(), 8.92893_f32);
    }

    #[test]
    fn eat_comment_at_eof_is_harmless() {
        let mut s = Scanner::from_text("5", "test");
        assert_eq!(s.next_i32("five").unwrap(), 5);
        s.eat_comment(); // C++ putback on a failed stream; here: no-op
        assert!(s.next_i32("past eof").is_err());
    }
}
