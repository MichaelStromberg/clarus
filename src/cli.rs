//! Shared CLI output helpers for Clarus binaries.

use std::time::Instant;

use colored::Colorize;

use crate::perf;

/// Print the Clarus banner with a subtitle line.
pub fn banner(subtitle: &str) {
    eprintln!();
    eprintln!("{} {}", "Clarus".bold().cyan(), subtitle.dimmed());
    eprintln!("{}", "(c) 2026 Michael Stromberg".dimmed());
    eprintln!();
}

/// Print a bold section header with a horizontal rule.
pub fn section(title: &str) {
    let bar = "─".repeat(50);
    eprintln!("{} {}", title.bold().blue(), bar.dimmed());
}

/// Print a key-value pair in aligned columns.
pub fn kv(key: &str, value: &str) {
    eprintln!("  {:<22} {}", key.dimmed(), value);
}

/// Format a number with thousands separators (e.g., 207367 → "207,367").
pub fn num(n: impl std::fmt::Display) -> String {
    let s = n.to_string();
    let mut result = String::with_capacity(s.len() + s.len() / 3);
    for (i, c) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 && c != '-' {
            result.push(',');
        }
        result.push(c);
    }
    result.chars().rev().collect()
}

/// Print a success message with a green checkmark.
pub fn success(msg: &str) {
    eprintln!("  {} {}", "✓".green().bold(), msg);
}

/// Print a warning message in yellow.
pub fn warning(msg: &str) {
    eprintln!("  {} {}", "⚠".yellow(), msg.yellow());
}

/// Print elapsed time and peak memory usage.
pub fn print_summary(start: Instant) {
    let elapsed = start.elapsed();
    eprintln!();
    eprintln!(
        "{}  {}\n{}  {}",
        "Time".dimmed(),
        perf::format_elapsed(elapsed).bold(),
        "Peak memory".dimmed(),
        perf::peak_memory_bytes()
            .map_or_else(|| "N/A".to_string(), perf::format_bytes)
            .bold(),
    );
    eprintln!();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_num() {
        assert_eq!(num(0), "0");
        assert_eq!(num(999), "999");
        assert_eq!(num(1000), "1,000");
        assert_eq!(num(1234567), "1,234,567");
    }
}
