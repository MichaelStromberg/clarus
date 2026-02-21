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

    Ok(filename.to_string())
}

/// Return `/tmp/{filename}` for a given URL.
pub fn local_path_for_url(url: &str) -> Result<PathBuf> {
    let filename = filename_from_url(url)?;
    Ok(PathBuf::from("/tmp").join(filename))
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
pub fn download_file(url: &str, dest: &Path) -> Result<()> {
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
    pub expected_md5: String,
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

impl ChecksumFailure {
    pub fn name(&self) -> &str {
        match self {
            Self::Mismatch { name, .. } | Self::ReadError { name, .. } => name,
        }
    }
}

/// Download multiple files in parallel using std::thread.
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

/// Verify MD5 checksums for a batch of files. Returns all failures.
pub fn verify_checksums(tasks: &[DownloadTask]) -> Vec<ChecksumFailure> {
    tasks
        .iter()
        .filter_map(|task| {
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
            if actual != task.expected_md5 {
                Some(ChecksumFailure::Mismatch {
                    name: task.name.clone(),
                    path: task.dest.clone(),
                    expected: task.expected_md5.clone(),
                    actual,
                })
            } else {
                None
            }
        })
        .collect()
}

/// Resolve files: check existence, verify MD5, delete on mismatch, return tasks to download.
///
/// Uses `compute_md5` directly instead of a separate existence check to avoid TOCTOU races.
/// Stale file removal is best-effort since `download_file` uses `File::create` which truncates.
pub fn resolve_files(entries: &[(&str, &str, &str)]) -> Result<Vec<DownloadTask>> {
    let mut to_download = Vec::new();

    for &(name, url, expected_md5) in entries {
        let local = local_path_for_url(url)?;

        let needs_download = match compute_md5(&local) {
            Ok(actual) if actual == expected_md5 => false,
            Ok(_) => {
                // MD5 mismatch — attempt to remove stale file (non-fatal if removal fails,
                // since download_file uses File::create which truncates)
                let _ = fs::remove_file(&local);
                true
            }
            Err(_) => true, // File missing or unreadable
        };

        if needs_download {
            to_download.push(DownloadTask {
                name: name.to_string(),
                url: url.to_string(),
                dest: local,
                expected_md5: expected_md5.to_string(),
            });
        }
    }

    Ok(to_download)
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
    fn local_path_for_url_returns_tmp() {
        let path = local_path_for_url("https://example.com/genome.fna.gz").unwrap();
        assert_eq!(path, PathBuf::from("/tmp/genome.fna.gz"));
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
            expected_md5: "eb733a00c0c9d336e65691a37ab54293".to_string(),
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
            expected_md5: "0000000000000000000000000000dead".to_string(),
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
            expected_md5: "0000000000000000000000000000dead".to_string(),
        }];

        let failures = verify_checksums(&tasks);
        assert_eq!(failures.len(), 1);
        assert!(
            matches!(&failures[0], ChecksumFailure::ReadError { name, .. } if name == "missing")
        );
    }
}
