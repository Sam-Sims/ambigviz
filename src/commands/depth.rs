use rust_htslib::bam::{IndexedReader, Read};

use crate::commands::utils::{get_start_stop, get_tids, read_bam};

pub struct Depth<'a> {
    input: &'a str,
    chrom: Option<&'a str>,
    start: u32,
    stop: u32,
    output: String,
}

impl<'a> Depth<'a> {
    pub fn new(
        input: &'a str,
        chrom: Option<&'a str>,
        start: Option<u32>,
        stop: Option<u32>,
        output: String,
    ) -> Self {
        let (start, stop) = get_start_stop(start, stop);
        Self {
            input,
            chrom,
            start,
            stop,
            output,
        }
    }

    fn plot(&self, x: Vec<u32>, y: Vec<u32>) {
        let trace = plotly::Scatter::new(x, y)
            .name("Depth")
            .mode(plotly::common::Mode::Lines);

        let mut plot = plotly::Plot::new();
        plot.add_trace(trace);
        let out_name = format!("{}.png", self.output);
        plot.write_image(out_name, plotly::ImageFormat::PNG, 2000, 1000, 1.0);
    }

    fn process_pileup(&self, bam: &mut IndexedReader) -> (Vec<u32>, Vec<u32>) {
        let mut x = Vec::new();
        let mut y = Vec::new();

        for p in bam.pileup().flatten() {
            let pos = p.pos();
            if pos >= self.start && pos < self.stop {
                x.push(pos + 1);
                y.push(p.depth());
            }
        }
        (x, y)
    }

    pub fn run(&self) {
        let mut bam = read_bam(self.input);
        let tids = get_tids(self.chrom, Some(bam.header()));
        println!("Tids: {:?}", tids);

        for tid in &tids {
            bam.fetch((tid, self.start, self.stop))
                .expect("Failed to fetch region");
            let (x, y) = self.process_pileup(&mut bam);
            self.plot(x, y);
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use rust_htslib::bam::{self, Record};

    use super::*;

    #[fixture]
    fn testbam() -> IndexedReader {
        let temp_path = "testing/depth_test.bam";

        // first we build a header
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
            Record::from_sam(
                &header_view,
                b"read1\t3\tchr1\t5\t60\t9M1D\tchr1\t80\t10\tGGGGGGGGG\tFFFFFFFFF",
            )
            .unwrap(),
            // Insertion of AA at end
            Record::from_sam(
                &header_view,
                b"read2\t3\tchr1\t5\t60\t10M2I\tchr1\t80\t10\tGGGGGGGGGGAA\tFFFFFFFFFFFF",
            )
            .unwrap(),
            // Insertion of AA at end
            Record::from_sam(
                &header_view,
                b"read3\t3\tchr1\t5\t60\t10M2I\tchr1\t80\t10\tGGGGGGGGGGAA\tFFFFFFFFFFFF",
            )
            .unwrap(),
            Record::from_sam(
                &header_view,
                b"read4\t3\tchr1\t5\t60\t10M\tchr1\t80\t10\tAGGGGGGGGG\tFFFFFFFFFF",
            )
            .unwrap(),
            Record::from_sam(
                &header_view,
                b"read5\t3\tchr1\t5\t60\t10M\tchr1\t80\t10\tAGGGGGGGGG\tFFFFFFFFFF",
            )
            .unwrap(),
            Record::from_sam(
                &header_view,
                b"read6\t3\tchr1\t5\t60\t7M\tchr1\t80\t10\tAGGGGGG\tFFFFFFF",
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
    fn test_process_pileup() {
        let mut bam = testbam();
        let depth = Depth::new(
            "testing/depth_test.bam",
            Some("chr1"),
            Some(1),
            Some(15),
            "testing/depth_test".to_string(),
        );
        let (x, y) = depth.process_pileup(&mut bam);
        assert_eq!(x, vec![5, 6, 7, 8, 9, 10, 11, 12, 13, 14]);
        assert_eq!(y, vec![6, 6, 6, 6, 6, 6, 6, 5, 5, 5]);
    }
}
