//! Shared CLI output helpers for Clarus binaries.

use std::time::Instant;

use colored::Colorize;

use crate::perf;

pub fn banner(subtitle: &str) {
    eprintln!();
    eprintln!("{} {}", "Clarus".bold().cyan(), subtitle.dimmed());
    eprintln!("{}", "(c) 2026 Michael Stromberg".dimmed());
    eprintln!();
}

pub fn section(title: &str) {
    let bar = "─".repeat(50);
    eprintln!("{} {}", title.bold().blue(), bar.dimmed());
}

pub fn kv(key: &str, value: &str) {
    eprintln!("  {:<20} {}", key.dimmed(), value);
}

pub fn success(msg: &str) {
    eprintln!("  {} {}", "✓".green().bold(), msg);
}

pub fn warning(msg: &str) {
    eprintln!("  {} {}", "⚠".yellow(), msg.yellow());
}

pub fn print_summary(start: Instant) {
    let elapsed = start.elapsed();
    eprintln!();
    eprintln!(
        "{}  {}\n{}  {}",
        "Time".dimmed(),
        perf::format_elapsed(elapsed).bold(),
        "Peak memory".dimmed(),
        perf::peak_memory_bytes()
            .map(perf::format_bytes)
            .unwrap_or_else(|| "N/A".to_string())
            .bold(),
    );
    eprintln!();
}
