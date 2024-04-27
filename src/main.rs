pub use crate::cli::Cli;
use clap::Parser;
mod cli;
pub mod commands;

fn main() {
    let args = Cli::parse();
    match args.command {
        cli::Command::Depth {
            input,
            region,
            output,
        } => {
            if let Some(region) = region {
                let region_string = cli::parse_region(&region);
                let (chrom, start, stop) = region_string;
                let depth_plotter =
                    commands::depth::Depth::new(&input, Some(chrom), Some(start), stop, output);
                depth_plotter.run();
            } else {
                let depth_plotter = commands::depth::Depth::new(&input, None, None, None, output);
                depth_plotter.run();
            }
        }
        cli::Command::Ambig {
            input,
            region,
            no_indel,
            threshold,
            no_label,
            output,
            base_quality_threshold,
            map_quality_threshold,
            depth_threshold,
            minor_depth_threshold,
            strand_bias_threshold,
            bed
        } => {
            if let Some(region) = region {
                let region_string = cli::parse_region(&region);
                let (chrom, start, stop) = region_string;
                let _plotter = commands::ambig::Ambig::new(
                    &input,
                    Some(chrom),
                    Some(start),
                    stop,
                    no_indel,
                    threshold,
                    no_label,
                    output,
                    base_quality_threshold,
                    map_quality_threshold,
                    depth_threshold,
                    minor_depth_threshold,
                    strand_bias_threshold,
                    bed
                );
                _plotter.run();
            } else {
                println!("No region specified");
                let ambig_plotter = commands::ambig::Ambig::new(
                    &input,
                    None,
                    None,
                    None,
                    no_indel,
                    threshold,
                    no_label,
                    output,
                    base_quality_threshold,
                    map_quality_threshold,
                    depth_threshold,
                    minor_depth_threshold,
                    strand_bias_threshold,
                    bed
                );
                ambig_plotter.run();
            }
        }
    }
}

//todo
// 1. Output nice tsv
// 2. Tests
