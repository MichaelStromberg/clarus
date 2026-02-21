use std::fs::{self, File};
use std::io::{BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::thread;

use anyhow::{Context, Result, bail};
use md5::{Digest, Md5};

const BUFFER_SIZE: usize = 64 * 1024;

/// Extract the filename from a URL's last path segment.
pub fn filename_from_url(url: &str) -> Result<String> {
    let path = url
        .split('?')
        .next()
        .unwrap_or(url)
        .split('#')
        .next()
        .unwrap_or(url);

    let filename = path.rsplit('/').next().unwrap_or("");

    if filename.is_empty() {
        bail!("cannot extract filename from URL: '{url}'");
    }

    if filename.contains('/') || filename.contains('\\') || filename == ".." {
        bail!("unsafe filename extracted from URL: '{filename}'");
    }

    Ok(filename.to_string())
}

/// Compute the MD5 hex digest of a file by streaming 64KB chunks.
pub fn compute_md5(path: &Path) -> Result<String> {
    let file = File::open(path)
        .with_context(|| format!("failed to open file for MD5: {}", path.display()))?;
    let mut reader = BufReader::new(file);
    let mut hasher = Md5::new();
    let mut buffer = vec![0u8; BUFFER_SIZE];

    loop {
        let n = reader.read(&mut buffer)?;
        if n == 0 {
            break;
        }
        hasher.update(&buffer[..n]);
    }

    Ok(format!("{:x}", hasher.finalize()))
}

/// Download a file from a URL via HTTP GET, streaming to disk.
fn download_file(url: &str, dest: &Path) -> Result<()> {
    let mut response = ureq::get(url)
        .call()
        .with_context(|| format!("HTTP request failed for {url}"))?;

    let mut file =
        File::create(dest).with_context(|| format!("failed to create {}", dest.display()))?;
    let mut buffer = vec![0u8; BUFFER_SIZE];
    let mut reader = response.body_mut().as_reader();

    loop {
        let n = reader.read(&mut buffer)?;
        if n == 0 {
            break;
        }
        file.write_all(&buffer[..n])?;
    }

    file.sync_all()?;
    Ok(())
}

/// A file to download.
#[derive(Debug, Clone)]
pub struct DownloadTask {
    pub name: String,
    pub url: String,
    pub dest: PathBuf,
    /// `Some(hash)` — verify MD5 after download; `None` — no MD5 available, existence-only caching.
    pub expected_md5: Option<String>,
}

/// Result of a single download attempt.
#[derive(Debug)]
pub struct DownloadResult {
    pub name: String,
    pub error: Option<String>,
}

/// A checksum verification failure.
#[derive(Debug)]
pub enum ChecksumFailure {
    Mismatch {
        name: String,
        path: PathBuf,
        expected: String,
        actual: String,
    },
    ReadError {
        name: String,
        path: PathBuf,
        error: String,
    },
}

/// Download multiple files in parallel using std::thread.
#[must_use]
pub fn download_parallel(tasks: &[DownloadTask]) -> Vec<DownloadResult> {
    let handles: Vec<_> = tasks
        .iter()
        .map(|task| {
            let name = task.name.clone();
            let url = task.url.clone();
            let dest = task.dest.clone();
            let handle = thread::spawn(move || {
                let error = match download_file(&url, &dest) {
                    Ok(()) => None,
                    Err(e) => Some(format!("{e:#}")),
                };
                DownloadResult { name, error }
            });
            (task.name.clone(), handle)
        })
        .collect();

    handles
        .into_iter()
        .map(|(name, handle)| match handle.join() {
            Ok(result) => result,
            Err(payload) => {
                let msg = payload
                    .downcast_ref::<&str>()
                    .map(|s| s.to_string())
                    .or_else(|| payload.downcast_ref::<String>().cloned())
                    .unwrap_or_else(|| "unknown panic".to_string());
                DownloadResult {
                    name,
                    error: Some(format!("thread panicked: {msg}")),
                }
            }
        })
        .collect()
}

/// Verify MD5 checksums for a batch of files. Skips tasks with no expected MD5. Returns all failures.
#[must_use]
pub fn verify_checksums(tasks: &[DownloadTask]) -> Vec<ChecksumFailure> {
    tasks
        .iter()
        .filter(|task| task.expected_md5.is_some())
        .filter_map(|task| {
            let expected = task.expected_md5.as_ref().unwrap();
            let actual = match compute_md5(&task.dest) {
                Ok(md5) => md5,
                Err(e) => {
                    return Some(ChecksumFailure::ReadError {
                        name: task.name.clone(),
                        path: task.dest.clone(),
                        error: format!("{e:#}"),
                    });
                }
            };
            if actual != *expected {
                Some(ChecksumFailure::Mismatch {
                    name: task.name.clone(),
                    path: task.dest.clone(),
                    expected: expected.clone(),
                    actual,
                })
            } else {
                None
            }
        })
        .collect()
}

/// Resolve pre-built download tasks: check existence and MD5, return tasks that need downloading.
///
/// - `expected_md5: Some(hash)` — skip only if file exists AND MD5 matches; remove stale on mismatch
/// - `expected_md5: None` — skip if file exists (no MD5 available)
#[must_use]
pub fn resolve_tasks(tasks: &[DownloadTask]) -> Vec<DownloadTask> {
    let mut to_download = Vec::new();

    for task in tasks {
        let needs_download = match &task.expected_md5 {
            None => !task.dest.exists(),
            Some(expected) => match compute_md5(&task.dest) {
                Ok(actual) if actual == *expected => false,
                Ok(_) => {
                    let _ = fs::remove_file(&task.dest);
                    true
                }
                Err(_) => true,
            },
        };

        if needs_download {
            to_download.push(task.clone());
        }
    }

    to_download
}

/// Download files and verify checksums, reporting progress to stderr.
///
/// Creates destination directories, downloads in parallel, reports errors,
/// then verifies MD5 checksums. Returns `Ok(())` if all downloads and
/// checksums pass, or bails with an error description.
pub fn download_and_verify(tasks: &[DownloadTask]) -> Result<()> {
    use crate::cli;
    use colored::Colorize;

    // Create destination directories
    for task in tasks {
        if let Some(parent) = task.dest.parent() {
            fs::create_dir_all(parent)?;
        }
    }

    for task in tasks {
        cli::kv(&task.name, &task.url);
    }

    let results = download_parallel(tasks);
    let errors: Vec<_> = results.iter().filter(|r| r.error.is_some()).collect();

    if !errors.is_empty() {
        for err in &errors {
            eprintln!(
                "  {} {}: {}",
                "✗".red().bold(),
                err.name,
                err.error.as_ref().unwrap()
            );
        }
        bail!("{} download(s) failed", errors.len());
    }

    for result in &results {
        cli::success(&format!("{} downloaded", result.name));
    }

    let failures = verify_checksums(tasks);
    if !failures.is_empty() {
        for f in &failures {
            match f {
                ChecksumFailure::Mismatch {
                    name,
                    expected,
                    actual,
                    ..
                } => {
                    eprintln!(
                        "  {} {name}: expected {expected}, got {actual}",
                        "✗".red().bold()
                    );
                }
                ChecksumFailure::ReadError { name, error, .. } => {
                    eprintln!("  {} {name}: {error}", "✗".red().bold());
                }
            }
        }
        bail!("{} checksum verification(s) failed", failures.len());
    }

    cli::success("all checksums verified");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn filename_from_url_simple() {
        assert_eq!(
            filename_from_url("https://example.com/path/to/file.txt").unwrap(),
            "file.txt"
        );
    }

    #[test]
    fn filename_from_url_with_query() {
        assert_eq!(
            filename_from_url("https://example.com/path/file.gz?token=abc").unwrap(),
            "file.gz"
        );
    }

    #[test]
    fn filename_from_url_trailing_slash_fails() {
        assert!(filename_from_url("https://example.com/path/").is_err());
    }

    #[test]
    fn filename_from_url_rejects_path_traversal() {
        assert!(filename_from_url("https://example.com/..").is_err());
    }

    #[test]
    fn compute_md5_known_value() {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(b"hello world").unwrap();
        f.flush().unwrap();
        let md5 = compute_md5(f.path()).unwrap();
        assert_eq!(md5, "5eb63bbbe01eeed093cb22bb8f5acdc3");
    }

    #[test]
    fn verify_checksums_pass() {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(b"test data").unwrap();
        f.flush().unwrap();

        let tasks = vec![DownloadTask {
            name: "test".to_string(),
            url: "https://example.com/test".to_string(),
            dest: f.path().to_path_buf(),
            expected_md5: Some("eb733a00c0c9d336e65691a37ab54293".to_string()),
        }];

        let failures = verify_checksums(&tasks);
        assert!(failures.is_empty());
    }

    #[test]
    fn verify_checksums_mismatch() {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(b"test data").unwrap();
        f.flush().unwrap();

        let tasks = vec![DownloadTask {
            name: "test".to_string(),
            url: "https://example.com/test".to_string(),
            dest: f.path().to_path_buf(),
            expected_md5: Some("0000000000000000000000000000dead".to_string()),
        }];

        let failures = verify_checksums(&tasks);
        assert_eq!(failures.len(), 1);
        assert!(matches!(
            &failures[0],
            ChecksumFailure::Mismatch { name, .. } if name == "test"
        ));
    }

    #[test]
    fn verify_checksums_read_error() {
        let tasks = vec![DownloadTask {
            name: "missing".to_string(),
            url: "https://example.com/missing".to_string(),
            dest: PathBuf::from("/tmp/nonexistent_clarus_test_file"),
            expected_md5: Some("0000000000000000000000000000dead".to_string()),
        }];

        let failures = verify_checksums(&tasks);
        assert_eq!(failures.len(), 1);
        assert!(
            matches!(&failures[0], ChecksumFailure::ReadError { name, .. } if name == "missing")
        );
    }

    #[test]
    fn verify_checksums_empty_md5_reports_mismatch() {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(b"test data").unwrap();
        f.flush().unwrap();

        let tasks = vec![DownloadTask {
            name: "empty_md5".to_string(),
            url: "https://example.com/empty_md5".to_string(),
            dest: f.path().to_path_buf(),
            expected_md5: Some(String::new()),
        }];

        let failures = verify_checksums(&tasks);
        assert_eq!(failures.len(), 1);
        assert!(matches!(
            &failures[0],
            ChecksumFailure::Mismatch { expected, actual, .. }
                if expected.is_empty() && actual == "eb733a00c0c9d336e65691a37ab54293"
        ));
    }

    #[test]
    fn verify_checksums_skips_none_md5() {
        let tasks = vec![DownloadTask {
            name: "no_md5".to_string(),
            url: "https://example.com/no_md5".to_string(),
            dest: PathBuf::from("/tmp/nonexistent_clarus_test_file"),
            expected_md5: None,
        }];

        let failures = verify_checksums(&tasks);
        assert!(failures.is_empty());
    }

    #[test]
    fn resolve_tasks_missing_file() {
        let tasks = vec![DownloadTask {
            name: "missing".to_string(),
            url: "https://example.com/missing".to_string(),
            dest: PathBuf::from("/tmp/nonexistent_clarus_resolve_test"),
            expected_md5: Some("dba0915f2560e6b9d4943fac274816b9".to_string()),
        }];

        let to_download = resolve_tasks(&tasks);
        assert_eq!(to_download.len(), 1);
        assert_eq!(to_download[0].name, "missing");
    }

    #[test]
    fn resolve_tasks_cached_file() {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(b"test data").unwrap();
        f.flush().unwrap();

        let tasks = vec![DownloadTask {
            name: "cached".to_string(),
            url: "https://example.com/cached".to_string(),
            dest: f.path().to_path_buf(),
            expected_md5: Some("eb733a00c0c9d336e65691a37ab54293".to_string()),
        }];

        let to_download = resolve_tasks(&tasks);
        assert!(to_download.is_empty());
    }

    #[test]
    fn resolve_tasks_md5_mismatch() {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(b"stale data").unwrap();
        f.flush().unwrap();

        let tasks = vec![DownloadTask {
            name: "stale".to_string(),
            url: "https://example.com/stale".to_string(),
            dest: f.path().to_path_buf(),
            expected_md5: Some("0000000000000000000000000000dead".to_string()),
        }];

        let to_download = resolve_tasks(&tasks);
        assert_eq!(to_download.len(), 1);
        assert_eq!(to_download[0].name, "stale");
    }

    #[test]
    fn resolve_tasks_empty_md5_always_redownloads() {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(b"some content").unwrap();
        f.flush().unwrap();

        let tasks = vec![DownloadTask {
            name: "empty_md5".to_string(),
            url: "https://example.com/empty_md5".to_string(),
            dest: f.path().to_path_buf(),
            expected_md5: Some(String::new()),
        }];

        let to_download = resolve_tasks(&tasks);
        assert_eq!(to_download.len(), 1);
        assert_eq!(to_download[0].name, "empty_md5");
    }

    #[test]
    fn resolve_tasks_none_md5_cached() {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(b"some content").unwrap();
        f.flush().unwrap();

        let tasks = vec![DownloadTask {
            name: "no_md5".to_string(),
            url: "https://example.com/no_md5".to_string(),
            dest: f.path().to_path_buf(),
            expected_md5: None,
        }];

        let to_download = resolve_tasks(&tasks);
        assert!(to_download.is_empty());
    }

    #[test]
    fn resolve_tasks_none_md5_missing() {
        let tasks = vec![DownloadTask {
            name: "no_md5_missing".to_string(),
            url: "https://example.com/no_md5_missing".to_string(),
            dest: PathBuf::from("/tmp/nonexistent_clarus_none_md5_test"),
            expected_md5: None,
        }];

        let to_download = resolve_tasks(&tasks);
        assert_eq!(to_download.len(), 1);
        assert_eq!(to_download[0].name, "no_md5_missing");
    }
}
