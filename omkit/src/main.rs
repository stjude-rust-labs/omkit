//! `omkit`—a cross-platform command line toolkit for omics-based analysis.

use anyhow::Result;
use clap::Parser;

fn main() -> Result<()> {
    Cli::parse();
    Ok(())
}

/// A cross-platform command line toolkit for omics-based analysis.
#[derive(Debug, Parser)]
#[command(name = "omkit", version, about)]
pub struct Cli {}
