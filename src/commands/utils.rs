use std::error::Error;
use std::path;

use rust_htslib::bam;
use rust_htslib::bam::IndexedReader;

pub fn get_tids(chrom: Option<&str>, header: Option<&bam::HeaderView>) -> Vec<String> {
    match chrom {
        Some(chrom) => vec![chrom.to_string()],
        None => header
            .unwrap()
            .target_names()
            .iter()
            .map(|x| String::from_utf8(x.to_vec()).unwrap())
            .collect(),
    }
}

pub fn get_start_stop(start: Option<u32>, stop: Option<u32>) -> (u32, u32) {
    let coord = 1;
    match (start, stop) {
        (Some(s), Some(e)) => (s - coord, e),
        (Some(s), None) => (s - coord, u32::MAX),
        (None, None) => (0, u32::MAX),
        _ => panic!("Invalid start/stop combination"),
    }
}

fn build_index(path: &str) -> Result<(), Box<dyn Error>> {
    match bam::index::build(path, None, bam::index::Type::Bai, 1) {
        Ok(_) => Ok(()),
        Err(e) => Err(e.into()),
    }
}

pub fn read_bam(path: &str) -> IndexedReader {
    let index_path = format!("{}.bai", path);
    if !path::Path::new(&index_path).exists() {
        println!("No index found for {}", path);
        println!("Attempting to build index");
        match build_index(path) {
            Ok(_) => {
                println!("Index built successfully.");
            }
            Err(err) => {
                eprintln!("Error building the index: {}", err);
            }
        }
    }
    let reader = bam::IndexedReader::from_path(path).expect("Failed to open BAM file");
    reader
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use rust_htslib::bam::Read;

    use super::*;

    #[fixture]
    fn testbam() -> bam::Reader {
        let temp_path = "testing/util_test.bam";
        let mut header = bam::Header::new();
        header.push_record(
            &bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"chr1")
                .push_tag(b"LN", &1000),
        );
        let header_view = bam::HeaderView::from_header(&header);
        let mut writer = bam::Writer::from_path(temp_path, &header, bam::Format::Bam).unwrap();
        let record = bam::Record::from_sam(
            &header_view,
            b"read1\t3\tchr1\t5\t60\t10M\tchr1\t80\t10\tGGGGGGGGGG\tFFFFFFFFFF",
        )
        .unwrap();
        writer.write(&record).unwrap();
        drop(writer);
        bam::Reader::from_path(temp_path).unwrap()
    }

    #[rstest]
    fn test_tids_from_bam(testbam: bam::Reader) {
        let tids = get_tids(None, Some(testbam.header()));
        let expected_tids = vec!["chr1".to_string()];
        assert_eq!(tids, expected_tids);
    }

    #[rstest]
    fn test_tids_from_supplied() {
        let tids = get_tids(Some("chr1"), None);
        let expected_tids = vec!["chr1".to_string()];
        assert_eq!(tids, expected_tids);
    }

    #[rstest]
    fn test_start_stop_both() {
        let start = Some(1);
        let stop = Some(1000);
        let (start, stop) = get_start_stop(start, stop);
        assert_eq!(start, 0);
        assert_eq!(stop, 1000);
    }

    #[rstest]
    fn test_start_stop_start() {
        let start = Some(1);
        let stop = None;
        let (start, stop) = get_start_stop(start, stop);
        let expected_stop = u32::MAX;
        assert_eq!(start, 0);
        assert_eq!(stop, expected_stop);
    }

    #[rstest]
    fn test_start_stop_none() {
        let start = None;
        let stop = None;
        let (start, stop) = get_start_stop(start, stop);
        let expected_stop = u32::MAX;
        assert_eq!(start, 0);
        assert_eq!(stop, expected_stop);
    }

    #[rstest]
    fn test_build_index_success(testbam: bam::Reader) {
        let path = "testing/util_test.bam";
        let result = build_index(path);
        assert!(result.is_ok());
    }

    #[rstest]
    fn test_read_bam(testbam: bam::Reader) {
        let bam = read_bam("testing/util_test.bam");
        let tid = String::from_utf8(bam.header().target_names()[0].to_vec()).unwrap();
        assert_eq!(tid, "chr1");
    }
}
