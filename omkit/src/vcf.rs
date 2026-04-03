//! The `vcf` subcommand.

use anyhow::Result;
use clap::Parser;
use clap::Subcommand;

/// Arguments for VCF subcommands.
#[derive(Debug, Parser)]
pub struct Args {
    /// The VCF subcommand to run.
    #[command(subcommand)]
    pub command: Command,
}

/// VCF subcommands.
#[derive(Debug, Subcommand)]
pub enum Command {
    /// Lifts over a VCF file from one genome build to another.
    Liftover,

    /// Normalizes a VCF file.
    Normalize,
}

/// Executes the `vcf` subcommand.
pub fn execute(args: Args) -> Result<()> {
    match args.command {
        Command::Liftover => todo!(),
        Command::Normalize => todo!(),
    }
}
