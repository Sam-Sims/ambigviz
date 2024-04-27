use std::collections::BTreeMap;
use std::fs::File;
use std::io::Write;

use plotly::{
    Bar,
    common::{Marker, TextAnchor, TextPosition, Title},
    ImageFormat, layout::{Axis, AxisType::Category, BarMode},
};
use rust_htslib::{bam, bam::Read};

use crate::commands::utils::{get_start_stop, get_tids, read_bam};

struct Pileup {
    // (total, (forward, reverse))
    a: (u32, (u32, u32)),
    t: (u32, (u32, u32)),
    c: (u32, (u32, u32)),
    g: (u32, (u32, u32)),
    del: (u32, (u32, u32)),
    ins: (u32, (u32, u32)),
}

impl Pileup {
    fn new() -> Self {
        Self {
            a: (0, (0, 0)),
            t: (0, (0, 0)),
            c: (0, (0, 0)),
            g: (0, (0, 0)),
            del: (0, (0, 0)),
            ins: (0, (0, 0)),
        }
    }

    fn get_total(&self) -> u32 {
        self.a.0 + self.t.0 + self.c.0 + self.g.0 + self.del.0 + self.ins.0
    }

    fn get_major_variant(&self) -> char {
        let (major_variant, _) = [
            ('A', self.a.0),
            ('T', self.t.0),
            ('C', self.c.0),
            ('G', self.g.0),
            ('-', self.del.0),
            ('+', self.ins.0),
        ]
            .iter()
            .cloned() // Clone the values to avoid borrowing issues
            .max_by_key(|&(_, i)| i)
            .expect("No bases found in base_counts.");
        major_variant
    }

    fn get_strand_ratio(&self, base: char) -> f64 {
        let (total, (forward, reverse)) = match base {
            'A' => self.a,
            'T' => self.t,
            'C' => self.c,
            'G' => self.g,
            '-' => self.del,
            '+' => self.ins,
            _ => panic!("Invalid base"),
        };
        if total == 0 {
            return -1.0;
        }
        let forward_ratio = forward as f64 / total as f64;
        //let reverse_ratio = reverse as f64 / total as f64;
        forward_ratio
    }

    fn get_base_count(&self, base: char) -> u32 {
        match base {
            'A' => self.a.0,
            'T' => self.t.0,
            'C' => self.c.0,
            'G' => self.g.0,
            '-' => self.del.0,
            '+' => self.ins.0,
            _ => panic!("Invalid base"),
        }
    }

    fn is_ambiguous(&self) -> bool {
        let mut count = 0;
        if self.a.0 > 0 {
            count += 1;
        }
        if self.t.0 > 0 {
            count += 1;
        }
        if self.c.0 > 0 {
            count += 1;
        }
        if self.g.0 > 0 {
            count += 1;
        }
        if self.del.0 > 0 {
            count += 1;
        }
        if self.ins.0 > 0 {
            count += 1;
        }
        if count > 1 {
            return true;
        }
        false
    }
}

pub struct Ambig<'a> {
    input: &'a str,
    chrom: Option<&'a str>,
    start: u32,
    stop: u32,
    no_indel: bool,
    threshold: f64,
    no_label: bool,
    output: String,
    base_quality_threshold: u8,
    map_quality_threshold: u8,
    depth_threshold: u32,
    minor_depth_treshold: u32,
    strand_bias_threshold: f64,
    bed: bool,
}

impl<'a> Ambig<'a> {
    pub fn new(
        input: &'a str,
        chrom: Option<&'a str>,
        start: Option<u32>,
        stop: Option<u32>,
        no_indel: bool,
        threshold: f64,
        no_label: bool,
        output: String,
        base_quality_threshold: u8,
        map_quality_threshold: u8,
        depth_threshold: u32,
        minor_depth_treshold: u32,
        strand_bias_threshold: f64,
        bed: bool
    ) -> Self {
        let (start, stop) = get_start_stop(start, stop);
        Self {
            input,
            chrom,
            start,
            stop,
            no_indel,
            threshold,
            no_label,
            output,
            base_quality_threshold,
            map_quality_threshold,
            depth_threshold,
            minor_depth_treshold,
            strand_bias_threshold,
            bed,
        }
    }

    fn output_tsv(&self, pos_to_plot: &BTreeMap<u32, BTreeMap<char, f64>>) {
        let mut file = File::create("output.tsv").expect("Failed to create file");
        // Iterate through the outer map
        for (key1, inner_map) in pos_to_plot {
            write!(file, "{}\t", key1).expect("Failed to write to file");
            // Iterate through the inner map
            for (key2, value) in inner_map {
                write!(file, "{}\t{}\t", key2, value).expect("Failed to write to file");
            }
            // Write a newline character to start a new row
            writeln!(file, "").expect("Failed to write to file");
        }
    }

    fn create_bar(
        &self,
        name: &str,
        colour: &str,
        pos: Vec<u32>,
        bases: Vec<f64>,
    ) -> Box<Bar<u32, f64>> {
        let cloned_colour = colour.to_string();
        let mut bar = Bar::new(pos, bases.clone()).name(name);
        if !self.no_label {
            bar = bar
                .text_array(
                    bases
                        .iter()
                        .map(|x| format!("{:.2}", x))
                        .collect::<Vec<String>>(),
                )
                .text_position(TextPosition::Inside)
                .inside_text_anchor(TextAnchor::Middle);
        }
        bar.marker(Marker::new().color(cloned_colour))
    }

    fn plot(&self, pos_to_plot: &BTreeMap<u32, BTreeMap<char, f64>>, tid: &str) {
        // collect all posisitons for x-axis
        let pos: Vec<u32> = pos_to_plot.keys().cloned().collect();

        let mut a: Vec<f64> = Vec::new();
        let mut c: Vec<f64> = Vec::new();
        let mut g: Vec<f64> = Vec::new();
        let mut t: Vec<f64> = Vec::new();
        let mut del: Vec<f64> = Vec::new();
        let mut ins: Vec<f64> = Vec::new();

        for base_counts in pos_to_plot.values() {
            a.push(*base_counts.get(&'A').unwrap_or(&0.0));
            c.push(*base_counts.get(&'C').unwrap_or(&0.0));
            g.push(*base_counts.get(&'G').unwrap_or(&0.0));
            t.push(*base_counts.get(&'T').unwrap_or(&0.0));
            del.push(*base_counts.get(&'-').unwrap_or(&0.0));
            ins.push(*base_counts.get(&'+').unwrap_or(&0.0));
        }

        let traces = vec![
            self.create_bar("A", "#60935D", pos.clone(), a.clone()),
            self.create_bar("C", "#1B5299", pos.clone(), c.clone()),
            self.create_bar("G", "#F5BB00", pos.clone(), g.clone()),
            self.create_bar("T", "#E63946", pos.clone(), t.clone()),
            self.create_bar("-", "#000000", pos.clone(), del.clone()),
            self.create_bar("+", "#6A041D", pos.clone(), ins.clone()),
        ];
        let layout = plotly::Layout::new()
            .bar_mode(BarMode::Stack)
            .title(Title::new("Ambiguous Bases"))
            .x_axis(Axis::new().title(Title::new("Position")).type_(Category))
            .y_axis(Axis::new().title(Title::new("Proportion")));
        let mut plot = plotly::Plot::new();
        for trace in traces {
            plot.add_trace(trace);
        }
        plot.set_layout(layout);

        let out_name = format!("{}_{}.png", tid, self.output);
        plot.write_image(out_name, ImageFormat::PNG, 2000, 1000, 1.0);
    }

    fn filter_base_counts(&self, pos: u32, pileup: &Pileup) -> BTreeMap<u32, BTreeMap<char, f64>> {
        let mut global_base_counts: BTreeMap<u32, BTreeMap<char, f64>> = BTreeMap::new();

        // First sum all the bases so we can calculate the percent later
        let total_count: u32 = pileup.get_total();

        // First find the major variant (base with most reads)
        let major_variant = pileup.get_major_variant();

        // check the sum of minor variants is greater than minor_depth_threshold
        if total_count - pileup.get_base_count(major_variant) < self.minor_depth_treshold {
            return global_base_counts;
        }
        // check the strand ratio for each base
        let strand_ratios: BTreeMap<char, f64> = [
            ('A', pileup.get_strand_ratio('A')),
            ('C', pileup.get_strand_ratio('C')),
            ('G', pileup.get_strand_ratio('G')),
            ('T', pileup.get_strand_ratio('T')),
            ('-', pileup.get_strand_ratio('-')),
            ('+', pileup.get_strand_ratio('+')),
        ]
            .iter()
            .cloned()
            .collect();
        for (base, ratio) in strand_ratios {
            // if ratio == -1.0 then there were no reads for that base
            if ratio == -1.0 {
                continue;
            }
            if base != major_variant && (ratio < self.strand_bias_threshold || ratio > 1.0-self.strand_bias_threshold){
                return global_base_counts;
            }
        }

        // Calculate the percent of each base and round to 4 decimal places, ignoring bases with 0 counts
        let percent_base_counts: BTreeMap<char, f64> = [
            ('A', pileup.a.0 as f64 / total_count as f64),
            ('C', pileup.c.0 as f64 / total_count as f64),
            ('G', pileup.g.0 as f64 / total_count as f64),
            ('T', pileup.t.0 as f64 / total_count as f64),
            ('-', pileup.del.0 as f64 / total_count as f64),
            ('+', pileup.ins.0 as f64 / total_count as f64),
        ]
            .iter()
            .cloned()
            .filter(|(_, percent)| *percent > 0.0)
            .map(|(base, percent)| (base, (percent * 10000.0).round() / 10000.0))
            .collect();

        // Find the proportion of minor variants
        let total_minor_proportion = percent_base_counts
            .iter()
            .filter(|(base, _)| *base != &major_variant)
            .map(|(_, percent)| percent)
            .sum::<f64>();

        // If the sum of the minor variant proportion is greater than the threshold, we will plot
        if total_minor_proportion > self.threshold {
            global_base_counts.insert(pos + 1, percent_base_counts);
        }
        global_base_counts
    }

    fn is_qc_pass(
        &self,
        record: &bam::Record,
        pileup: &bam::pileup::Pileup,
        alignment: &bam::pileup::Alignment,
    ) -> bool {
        // Checking the seq isnt empty is necessary as secondary alignments can cause empty seqs
        if record.seq().is_empty() || pileup.depth() < self.depth_threshold {
            return false;
        }
        // check for base Q score and map Q score
        if let Some(qpos) = alignment.qpos() {
            if record.qual()[qpos] < self.base_quality_threshold {
                return false;
            }
            if record.mapq() < self.map_quality_threshold {
                return false;
            }
        }
        true
    }

    fn produce_pileup(&self, bam: &mut bam::IndexedReader) -> BTreeMap<u32, BTreeMap<char, f64>> {
        let mut filtered_pileup_counts: BTreeMap<u32, BTreeMap<char, f64>> = BTreeMap::new();
        for pileup in bam.pileup().flatten() {
            let mut pileup_struct = Pileup::new();
            if pileup.pos() < self.start || pileup.pos() >= self.stop {
                continue;
            }
            //let mut base_counts: HashMap<char, u32> = HashMap::new();
            for alignment in pileup.alignments() {
                let record = alignment.record();
                let orientation = record.strand().to_string();
                if !self.is_qc_pass(&record, &pileup, &alignment) {
                    continue;
                }
                if self.no_indel && (alignment.is_refskip() || alignment.is_del()) {
                    continue;
                } else if alignment.is_del() {
                    pileup_struct.del.0 += 1;
                    if orientation == "+" {
                        pileup_struct.del.1.0 += 1;
                    } else {
                        pileup_struct.del.1.1 += 1;
                    }
                }
                // if read passes qc, is not a deletion or a refskip then we have a real base
                if let Some(qpos) = alignment.qpos() {
                    let base = (record.seq()[qpos] as char).to_ascii_uppercase();
                    match base {
                        'A' => {
                            pileup_struct.a.0 += 1;
                            if orientation == "+" {
                                pileup_struct.a.1.0 += 1;
                            } else {
                                pileup_struct.a.1.1 += 1;
                            }
                        }
                        'T' => {
                            pileup_struct.t.0 += 1;
                            if orientation == "+" {
                                pileup_struct.t.1.0 += 1;
                            } else {
                                pileup_struct.t.1.1 += 1;
                            }
                        }
                        'C' => {
                            pileup_struct.c.0 += 1;
                            if orientation == "+" {
                                pileup_struct.c.1.0 += 1;
                            } else {
                                pileup_struct.c.1.1 += 1;
                            }
                        }
                        'G' => {
                            pileup_struct.g.0 += 1;
                            if orientation == "+" {
                                pileup_struct.g.1.0 += 1;
                            } else {
                                pileup_struct.g.1.1 += 1;
                            }
                        }
                        _ => {}
                    }
                }
                // Check if insertion
                if let bam::pileup::Indel::Ins(_len) = alignment.indel() {
                    if self.no_indel {
                        continue;
                    } else {
                        pileup_struct.ins.0 += 1;
                        if orientation == "+" {
                            pileup_struct.ins.1.0 += 1;
                        } else {
                            pileup_struct.ins.1.1 += 1;
                        }
                    }
                }
            }

            // skip processing any if only 1 base present (no ambiguity) at that position
            if pileup_struct.is_ambiguous() {
                let pileup_counts = self.filter_base_counts(pileup.pos(), &pileup_struct);
                if !pileup_counts.is_empty() {
                    filtered_pileup_counts.extend(pileup_counts);
                }
            }
        }
        filtered_pileup_counts
    }

    fn output_bed(&self, tid: &String, pos_to_plot: BTreeMap<u32, BTreeMap<char, f64>>) {
        let out_name = format!("{}_{}.bed", tid, self.output);
        let mut file = File::create(out_name).expect("Failed to create file");
        for (pos, _) in pos_to_plot {
            writeln!(file, "{}\t{}\t{}\t{}", tid, pos, pos + 1, tid).expect("Failed to write to file");
        }
    }

    pub fn run(&self) {
        let mut bam = read_bam(self.input);
        let tids = get_tids(self.chrom, Some(bam.header()));
        println!("Tids: {:?}", tids);
        // run for each chromosome
        for tid in &tids {
            println!("Processing Tid: {}", tid);
            bam.fetch((tid, self.start, self.stop))
                .expect("Failed to fetch region");
            let filtered_pileup_counts = self.produce_pileup(&mut bam);
            self.plot(&filtered_pileup_counts, tid);
            //self.output_tsv(&pos_to_plot);
            if self.bed {
                self.output_bed(tid, filtered_pileup_counts)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use rust_htslib::bam::IndexedReader;

    use super::*;

    #[fixture]
    fn testbam() -> IndexedReader {
        let temp_path = "testing/ambig_test.bam";
        let mut header = bam::Header::new();
        header.push_record(
            &bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"chr1")
                .push_tag(b"LN", &1000),
        );
        let header_view = bam::HeaderView::from_header(&header);
        let mut writer = bam::Writer::from_path(temp_path, &header, bam::Format::Bam).unwrap();
        let records = vec![
            // Deletion at end
            bam::Record::from_sam(
                &header_view,
                b"read1\t3\tchr1\t5\t60\t9M1D\tchr1\t80\t10\tGGGGGGGGG\tFFFFFFFFF",
            )
                .unwrap(),
            // Insertion of AA at end
            bam::Record::from_sam(
                &header_view,
                b"read2\t3\tchr1\t5\t60\t10M2I\tchr1\t80\t10\tGGGGGGGGGGAA\tFFFFFFFFFFFF",
            )
                .unwrap(),
            // Insertion of AA at end
            bam::Record::from_sam(
                &header_view,
                b"read3\t3\tchr1\t5\t60\t10M2I\tchr1\t80\t10\tGGGGGGGGGGAA\tFFFFFFFFFFFF",
            )
                .unwrap(),
            bam::Record::from_sam(
                &header_view,
                b"read4\t3\tchr1\t5\t60\t10M\tchr1\t80\t10\tAGGGGGGGGG\tFFFFFFFFFF",
            )
                .unwrap(),
            bam::Record::from_sam(
                &header_view,
                b"read5\t3\tchr1\t5\t60\t10M\tchr1\t80\t10\tAGGGGGGGGG\tFFFFFFFFFF",
            )
                .unwrap(),
            bam::Record::from_sam(
                &header_view,
                b"read6\t3\tchr1\t5\t60\t10M\tchr1\t80\t10\tAGGGGGGGGG\tFFFFFFFFFF",
            )
                .unwrap(),
        ];
        for record in records {
            writer.write(&record).unwrap();
        }
        drop(writer);
        bam::index::build(temp_path, None, bam::index::Type::Bai, 1).unwrap();

        bam::IndexedReader::from_path(temp_path).unwrap()
    }

    #[rstest]
    fn test_pileup_below_threshold(testbam: IndexedReader) {
        let mut bam = testbam;
        let ambig = Ambig::new(
            "",
            Some("chr1"),
            Some(4),
            Some(6),
            true,
            0.5,
            false,
            "".to_string(),
            1,
            1,
            1,
            1,
            0.0,
            false,
        );
        let pos_to_plot = ambig.produce_pileup(&mut bam);
        let expected_pos = BTreeMap::new();
        assert_eq!(pos_to_plot, expected_pos);
    }

    #[rstest]
    fn test_pileup_above_threshold(testbam: IndexedReader) {
        let mut bam = testbam;
        let ambig = Ambig::new(
            "",
            Some("chr1"),
            Some(4),
            Some(6),
            true,
            0.2,
            false,
            "".to_string(),
            1,
            1,
            1,
            1,
            0.0,
            false,
        );
        let pos_to_plot = ambig.produce_pileup(&mut bam);

        let expected_pos = {
            let mut expected_bases = BTreeMap::new();
            expected_bases.insert('A', 0.5);
            expected_bases.insert('G', 0.5);
            let mut expected_pos = BTreeMap::new();
            expected_pos.insert(5, expected_bases);
            expected_pos
        };
        assert_eq!(pos_to_plot, expected_pos);
    }

    #[rstest]
    fn test_pileup_indel(testbam: IndexedReader) {
        let mut bam = testbam;
        let ambig = Ambig::new(
            "",
            Some("chr1"),
            Some(13),
            Some(15),
            false,
            0.1,
            false,
            "".to_string(),
            1,
            1,
            1,
            1,
            0.0,
            false,
        );
        let pos_to_plot = ambig.produce_pileup(&mut bam);
        let expected_pos = {
            let mut expected_bases = BTreeMap::new();
            expected_bases.insert('-', 0.125);
            expected_bases.insert('+', 0.25);
            expected_bases.insert('G', 0.625);
            let mut expected_pos = BTreeMap::new();
            expected_pos.insert(14, expected_bases);
            expected_pos
        };
        assert_eq!(pos_to_plot, expected_pos);
    }

    #[rstest]
    fn test_pileup_no_indels(testbam: IndexedReader) {
        let mut bam = testbam;
        let ambig = Ambig::new(
            "",
            Some("chr1"),
            Some(13),
            Some(15),
            true,
            0.1,
            false,
            "".to_string(),
            1,
            1,
            1,
            1,
            0.0,
            false
        );
        let pos_to_plot = ambig.produce_pileup(&mut bam);
        let expected_pos = BTreeMap::new();
        assert_eq!(pos_to_plot, expected_pos);
    }

    // #[rstest]
    // fn test_process_pileup() {
    //     let base_counts = {
    //         let mut base_counts = HashMap::new();
    //         base_counts.insert('A', 3);
    //         base_counts.insert('C', 3);
    //         base_counts.insert('G', 4);
    //         base_counts.insert('T', 0);
    //         base_counts
    //     };
    //     let ambig = Ambig::new(
    //         "",
    //         Some("chr1"),
    //         Some(4),
    //         Some(6),
    //         true,
    //         0.1,
    //         false,
    //         "".to_string(),
    //         1,
    //         1,
    //         1,
    //         1,
    //     );
    //     let pos_to_plot = ambig.filter_base_counts(1, base_counts);
    //     let expected_pos = {
    //         let mut expected_bases = BTreeMap::new();
    //         expected_bases.insert('A', 0.3);
    //         expected_bases.insert('C', 0.3);
    //         expected_bases.insert('G', 0.4);
    //         expected_bases.insert('T', 0.0);
    //         let mut expected_pos = BTreeMap::new();
    //         expected_pos.insert(2, expected_bases);
    //         expected_pos
    //     };
    //     assert_eq!(pos_to_plot, expected_pos);
    // }
}
