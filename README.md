[![Release](https://github.com/Sam-Sims/ambigviz/actions/workflows/release.yaml/badge.svg)](https://github.com/Sam-Sims/Kractor/actions/workflows/release.yaml)
![GitHub release (with filter)](https://img.shields.io/github/v/release/sam-sims/ambigviz)
![crates.io](https://img.shields.io/crates/v/ambigviz
)

# Ambigviz

ambigviz is a tool for rapidly scanning and visualising ambiguous/mixed bases at given positions in a
BAM file. It was initially written to examine intrahost diversity / co-infection in SARS-CoV-2 samples however it can be used for any BAM file (and scales well thanks to Rust). It uses strict filtering options by default to avoid sequencing artifacts and contamination/sequencing errors, The idea is to rapidly produce plots that are "presentation-ready".

It provides a simple command line interface and requires at minimum a BAM file only.

![example](img/plot.png)

## Installation

### Binaries:

Precompiled binaries for Linux, MacOS and Windows are attached to the latest release [0.1.0](https://github.com/Sam-Sims/ambigviz/releases/tag/v0.1.0)

### Cargo:

Requires [cargo](https://www.rust-lang.org/tools/install)

```
cargo install ambigviz
```

### Build from source:

#### Install rust toolchain:

To install please refer to the rust documentation: [docs](https://www.rust-lang.org/tools/install)

#### Clone the repository:

```bash
git clone https://github.com/Sam-Sims/ambigviz
```

#### Build and add to path:

```bash
cd ambigviz
cargo build --release
export PATH=$PATH:$(pwd)/target/release
```

All executables will be in the directory ambigviz/target/release.

## Usage

### Basic usage:

```bash
ambigviz ambig <path_to_bam> <region> [options]
```

At minimum, you need to provide a BAM file path. Ambigviz will look for a BAM index file (`.bai`) in the same directory with the pattern `[input].bai`. If no index is found, Ambigviz will attempt to create one for you. If indexing fails, you can index the BAM file yourself using samtools.

If no region is specified, the entire BAM file will be scanned for ambiguous bases. With the default 20% ambiguity threshold, this provides a quick way to scan a BAM file for potential sequencing errors, contamination, or co-infection, and flag regions of interest for further investigation.

Regions should be specified in the samtools format: `chr:start-end`, with 1-based coordinate positions.

The `--bed` option allows you to output the identified ambiguous positions to a BED file.

You can also generate a plot of the read depth across a region using the depth command:

```bash
ambigviz depth <path_to_bam> <region> [options]
```

### Options:

### Output

`-o, --output <output>` | Default: `ambig`

The `--output` option allows you to specify the base filename for the output file(s). Since a BAM file can contain data for multiple chromosomes, the full output filename will be formatted as `<chromosome>_<output>` for each chromosome present.

If this option is not provided, the default base filename used is ambig. So the output files would be named `<chromosome>_ambig` for each chromosome.

#### Threshold

`-t, --threshold <threshold>` | Default: `0.2`

The `--threshold` option sets the minimum allowed proportion of ambiguous bases at a given position for that position to be included in the output plot. The value provided should be between 0 and 0.5.
The default value is 0.2, meaning that any position where at least 20% of the reads show an ambiguous/mixed base will be plotted.

For example, if a position has 75 reads showing 'A' and 25 reads showing 'T', since the ambiguous proportion (25/100 = 0.25) exceeds the 0.2 threshold, that position will be included in the plot.

#### Indels

`--no-indel` | Default: `False`

By default, indels (insertions/deletions) are included when counting the number of ambiguous bases at each position.

The `--no-indel` option allows you to exclude indels from being counted as ambiguous bases.
This can be useful when working with noisier sequencing data that tends to have more indel errors, such as ONT.

#### Depth

`-d, --depth <depth>` | Default: `100`

The `--depth` option sets the minimum total read depth (number of reads) required at a position for that position to be included in the output plot. The default value is 100, meaning positions with fewer than 100 reads mapped will be excluded.

Increasing this value filters out low-coverage positions, while decreasing it allows lower-depth positions to be plotted.


#### Base quality

`-q, --base-quality <base-quality>` | Default: `20`

The `--base-quality` option sets the minimum Phred-scaled base quality score required for an individual read base to be considered for analysis at a given position. The default value is 20.

Any read bases with quality scores below this threshold will be ignored when calculating whether the position meets the ambiguity threshold for plotting.

#### Mapping quality

`-Q, --mapping-quality <mapping-quality>` | Default: `60`

The --mapping-quality option sets the minimum Phred-scaled mapping quality score required for a read to be included in the analysis at a given position. The default value is 60.

Any reads with mapping quality scores below this threshold will be completely ignored when analysing that position.

#### Strand bias

`-s, --strand-bias <strand-bias>` | Default: `0.1`

The `--strand-bias` option sets the minimum threshold for the balance of reads coming from the forward and reverse strands, in order for that position to be included in the plot. The value provided should be between 0 and 0.5. The default value is 0.1, meaning that at least 10% of the total reads from the ambiguous base at a given position must come from each strand. This allows you to filter out ambiguous reads which have a significant imbalance between the forward and reverse strand reads.

For example, if there are 80 reads of 'A' and 20 reads of 'C' at a position, the 'C' is identified as the ambiguous base. With a value of 0.1, at least 2 of the 'C' reads (10% of 20) must come from each strand.

A value of 0.5 requires an equal number of reads from both strands. Setting the value to 0 disables this strand bias filter entirely, allowing all positions to be included in the plot regardless of strand balance.

#### Bed file

`--bed` | Default: `False`

The --bed option allows you to output the identified ambiguous positions to a BED file format, in addition to the main graphical plot output.

By default, this option is disabled, meaning no BED file will be generated.

#### Labels

`--no-labels` | Default: `False`

By default the proportion of each base for each position are included as text annotations. This option will exclude
them.

