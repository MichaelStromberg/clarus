//! Performance measurement utilities.

use std::time::Duration;

/// Formats a duration as HH:MM:SS.d (tenths of a second).
#[must_use]
pub fn format_elapsed(d: Duration) -> String {
    let total_secs = d.as_secs();
    let hours = total_secs / 3600;
    let minutes = (total_secs % 3600) / 60;
    let seconds = total_secs % 60;
    let tenths = d.subsec_millis() / 100;
    format!("{hours:02}:{minutes:02}:{seconds:02}.{tenths}")
}

/// Returns peak resident set size in bytes, or None if unavailable.
#[must_use]
pub fn peak_memory_bytes() -> Option<u64> {
    #[cfg(any(target_os = "macos", target_os = "linux"))]
    {
        use std::mem::MaybeUninit;
        let mut usage = MaybeUninit::<libc::rusage>::uninit();
        // SAFETY: `getrusage` with `RUSAGE_SELF` and a valid pointer to an uninitialized
        // `rusage` struct is well-defined. The pointer is valid because `MaybeUninit`
        // guarantees proper alignment and size.
        let ret = unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) };
        if ret == 0 {
            // SAFETY: `getrusage` returned 0 (success), which guarantees the `rusage`
            // struct has been fully initialized.
            let usage = unsafe { usage.assume_init() };
            let bytes = if cfg!(target_os = "macos") {
                // macOS: ru_maxrss is in bytes
                usage.ru_maxrss as u64
            } else {
                // Linux: ru_maxrss is in kilobytes
                usage.ru_maxrss as u64 * 1024
            };
            return Some(bytes);
        }
    }
    None
}

/// Formats a byte count as a human-readable string (B, KB, MB, GB).
#[must_use]
pub fn format_bytes(bytes: u64) -> String {
    const KB: u64 = 1024;
    const MB: u64 = 1024 * 1024;
    const GB: u64 = 1024 * 1024 * 1024;

    if bytes >= GB {
        format!("{:.1} GB", bytes as f64 / GB as f64)
    } else if bytes >= MB {
        format!("{:.1} MB", bytes as f64 / MB as f64)
    } else if bytes >= KB {
        format!("{:.1} KB", bytes as f64 / KB as f64)
    } else {
        format!("{bytes} B")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_elapsed_basic() {
        assert_eq!(format_elapsed(Duration::from_millis(4400)), "00:00:04.4");
        assert_eq!(format_elapsed(Duration::from_secs(61)), "00:01:01.0");
        assert_eq!(format_elapsed(Duration::from_secs(3661)), "01:01:01.0");
    }

    #[test]
    fn format_bytes_units() {
        assert_eq!(format_bytes(500), "500 B");
        assert_eq!(format_bytes(1536), "1.5 KB");
        assert_eq!(format_bytes(10 * 1024 * 1024), "10.0 MB");
        assert_eq!(format_bytes(3 * 1024 * 1024 * 1024), "3.0 GB");
    }

    #[test]
    fn peak_memory_returns_value() {
        // On supported platforms, we should get a value
        if cfg!(any(target_os = "macos", target_os = "linux")) {
            let mem = peak_memory_bytes();
            assert!(mem.is_some());
            assert!(mem.unwrap() > 0);
        }
    }
}
