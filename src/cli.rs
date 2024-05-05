use std::str::FromStr;

use clap::{Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(version, about = "Visulise ambigous bases in a bam file.")]
pub struct Cli {
    #[clap(subcommand)]
    pub command: Command,
}

#[derive(Debug, Subcommand)]
pub enum Command {
    /// Plot Depths
    Depth {
        // Path to input bam
        #[arg(required = true, value_parser(check_input_exists))]
        input: String,

        // SAMtools region string
        #[arg(required = false)]
        region: Option<String>,

        // Output
        #[arg(short = 'o', long = "output", default_value = "ambig")]
        output: String,
    },
    /// Plot Ambigous bases
    Ambig {
        // Path to input bam
        #[arg(required = true, value_parser(check_input_exists))]
        input: String,

        // SAMtools region string
        #[arg(required = false)]
        region: Option<String>,

        // Output
        #[arg(short = 'o', long = "output", default_value = "ambig")]
        output: String,

        // Threshold for ambigous bases
        #[arg(
            short = 't',
            long = "threshold",
            default_value = "0.1",
            value_parser(check_threshold_valid)
        )]
        threshold: f64,

        // Threshold for base quality
        #[arg(short = 'q', long = "--min-BQ", default_value = "20")]
        base_quality_threshold: u8,

        // Threshold for map quality
        #[arg(short = 'Q', long = "--min-MQ", default_value = "60")]
        map_quality_threshold: u8,

        // Threshold for total depth
        #[arg(short = 'd', long = "--depth", default_value = "100")]
        depth_threshold: u32,

        // Threshold for depth of minor allele
        // #[arg(short = 'd', long = "--minor-depth", default_value = "20")]
        // minor_depth_threshold: u32,

        // Threshold for strand bias
        #[arg(
            short = 's',
            long = "--strand-bias",
            default_value = "0.1",
            value_parser(check_threshold_valid)
        )]
        strand_bias_threshold: f64,

        // Do not include indels
        #[arg(long = "no-indel")]
        no_indel: bool,

        // Do not include labels
        #[arg(long = "no-label")]
        no_label: bool,

        // Output bed file instead of plot
        #[arg(long = "bed")]
        bed: bool,
    },
}

fn check_input_exists(s: &str) -> Result<String, String> {
    if std::path::Path::new(s).exists() {
        Ok(s.to_string())
    } else {
        Err(format!("File does not exist: {}", s))
    }
}

fn check_threshold_valid(s: &str) -> Result<f64, String> {
    let threshold = f64::from_str(s).unwrap();
    if threshold >= 0.0 && threshold <= 0.5 {
        Ok(threshold)
    } else {
        Err(format!("Threshold must be between 0 and 0.5"))
    }
}

pub fn parse_region(region: &str) -> (&str, u32, Option<u32>) {
    // Check if we have just chromosome or both chromosome and region
    if region.contains(':') {
        let region_parts: Vec<&str> = region.split(':').collect();

        let chrom = region_parts[0];
        let positions = region_parts[1];

        let start: u32;
        let end: Option<u32>;
        // Check if we have a range of positions or only one
        if positions.contains('-') {
            let position_parts: Vec<u32> = positions
                .split('-')
                .map(|x| u32::from_str(x).unwrap())
                .collect();
            start = position_parts[0];
            end = Option::from(position_parts[1]);
        } else {
            start = u32::from_str(positions).unwrap();
            end = None;
        }
        (chrom, start, end)
    } else {
        // We only have a chromosome
        (region, 1, None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_region_chrom_only() {
        let region = "chr1";
        let (chrom, start, end) = parse_region(region);

        assert_eq!(chrom, "chr1");
        assert_eq!(start, 1);
        assert_eq!(end, None);
    }

    #[test]
    fn test_parse_region_chrom_and_single_position() {
        let region = "chr2:100";
        let (chrom, start, end) = parse_region(region);

        assert_eq!(chrom, "chr2");
        assert_eq!(start, 100);
        assert_eq!(end, None);
    }

    #[test]
    fn test_parse_region_chrom_and_range() {
        let region = "chr3:200-300";
        let (chrom, start, end) = parse_region(region);

        assert_eq!(chrom, "chr3");
        assert_eq!(start, 200);
        assert_eq!(end, Some(300));
    }
}
