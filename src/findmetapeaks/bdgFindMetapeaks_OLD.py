import argparse
from pathlib import Path
import shutil
import subprocess
import os

try:
    from . import __version__
except ImportError:
    __version__ = "dev"

try:
    import pybedtools
except ImportError as e:
    raise ImportError(
        "pybedtools is required. Install it with: python -m pip install --user pybedtools"
    ) from e


def parse_bed_like_line(line: str, file_path: Path, line_num: int):
    """
    Parse one BED-like line.

    Returns:
        None if the line is blank or a header line.
        list[str] of fields if the line is valid BED-like content.

    Raises:
        ValueError if the line is malformed.
    """
    stripped = line.strip()

    if not stripped:
        return None
    if stripped.startswith("#") or stripped.startswith("track") or stripped.startswith("browser"):
        return None

    # Accept either tabs or spaces in input, then normalize to tabs on output.
    fields = stripped.split()

    if len(fields) < 3:
        raise ValueError(
            f"{file_path}: line {line_num} has fewer than 3 columns ({len(fields)})"
        )

    chrom, start_str, end_str = fields[0], fields[1], fields[2]

    if chrom == "":
        raise ValueError(f"{file_path}: line {line_num} has an empty chromosome field")

    try:
        start = int(start_str)
    except ValueError:
        raise ValueError(f"{file_path}: line {line_num} has non-integer start: {start_str}")

    try:
        end = int(end_str)
    except ValueError:
        raise ValueError(f"{file_path}: line {line_num} has non-integer end: {end_str}")

    if start < 0:
        raise ValueError(f"{file_path}: line {line_num} has negative start: {start}")
    if end <= start:
        raise ValueError(
            f"{file_path}: line {line_num} has end <= start: start={start}, end={end}"
        )

    return fields



def is_bed_like_file(file_path: Path, max_check_lines: int = 50):
    """
    Quick content-based test for whether a file looks BED-like.
    We do not rely on the filename extension.
    """
    if not file_path.is_file():
        return False, "not a regular file"

    checked_records = 0

    try:
        with file_path.open("r") as handle:
            for line_num, line in enumerate(handle, start=1):
                try:
                    fields = parse_bed_like_line(line, file_path, line_num)
                except ValueError as e:
                    return False, str(e)

                if fields is None:
                    continue

                checked_records += 1
                if checked_records >= max_check_lines:
                    break

    except UnicodeDecodeError:
        return False, "not plain text"
    except OSError as e:
        return False, f"could not read file: {e}"

    if checked_records == 0:
        return False, "no interval records found"

    return True, None



def find_bed_like_files(input_dir: Path):
    """
    Scan top-level files in the input directory.
    Accept any file whose content looks BED-like, regardless of suffix.
    """
    accepted_files = []
    skipped_files = []

    for path in sorted(input_dir.iterdir()):
        if not path.is_file():
            continue
        if path.name.startswith("."):
            continue

        ok, reason = is_bed_like_file(path)
        if ok:
            accepted_files.append(path)
        else:
            skipped_files.append((path, reason))

    return accepted_files, skipped_files



def concatenate_interval_files(interval_files, output_bed: Path) -> int:
    """
    Concatenate BED-like files into one BED-like file.
    Validates every data line and normalizes whitespace to tabs.
    Returns the number of interval records written.
    """
    records_written = 0

    with output_bed.open("w") as out_handle:
        for interval_file in interval_files:
            with interval_file.open("r") as in_handle:
                for line_num, line in enumerate(in_handle, start=1):
                    fields = parse_bed_like_line(line, interval_file, line_num)
                    if fields is None:
                        continue

                    out_handle.write("\t".join(fields) + "\n")
                    records_written += 1

    return records_written



def load_genome_chroms(genome_path: Path):
    """
    Read chromosome names from a standard 2-column genome file.
    """
    chroms = []

    with genome_path.open("r") as handle:
        for line_num, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                continue

            fields = stripped.split()
            if len(fields) < 2:
                raise ValueError(
                    f"{genome_path}: line {line_num} has fewer than 2 columns"
                )

            chrom = fields[0]
            size_str = fields[1]

            try:
                size = int(size_str)
            except ValueError:
                raise ValueError(
                    f"{genome_path}: line {line_num} has non-integer chromosome size: {size_str}"
                )

            if size <= 0:
                raise ValueError(
                    f"{genome_path}: line {line_num} has non-positive chromosome size: {size}"
                )

            chroms.append(chrom)

    if not chroms:
        raise ValueError(f"No chromosome entries found in genome file: {genome_path}")

    return chroms



def load_bed_chroms(bed_path: Path):
    """
    Read chromosome names from a BED-like file.
    """
    chroms = set()

    with bed_path.open("r") as handle:
        for line_num, line in enumerate(handle, start=1):
            fields = parse_bed_like_line(line, bed_path, line_num)
            if fields is None:
                continue
            chroms.add(fields[0])

    return chroms



def main():
    parser = argparse.ArgumentParser(prog="bdgFindMetapeaks")

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input directory containing BED-like interval files (e.g. BED, narrowPeak, broadPeak)"
    )
    parser.add_argument(
        "-o", "--output-dir",
        required=True,
        help="Main output directory"
    )
    parser.add_argument(
        "-G", "--genome-file",
        required=True,
        help="Genome file matching the genome build used for the input interval files"
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=5.0,
        help="Signal cutoff used to detect candidate enriched regions (default: 5.0)",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=200,
        help="Minimum peak length in bp (default: 200)",
    )
    parser.add_argument(
        "--max-gap",
        type=int,
        default=30,
        help="Maximum gap for merging adjacent enriched segments in bp (default: 30)",
    )

    args = parser.parse_args()

    if shutil.which("macs2") is None:
        raise RuntimeError("macs2 is not installed or not on PATH")
    if shutil.which("bedtools") is None:
        raise RuntimeError("bedtools is not installed or not on PATH")
    if shutil.which("sort") is None:
        raise RuntimeError("sort is not available on PATH")

    input_dir = Path(args.input)
    main_output_dir = Path(args.output_dir)
    genome_path = Path(args.genome_file)

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    if not input_dir.is_dir():
        raise ValueError(f"Input path is not a directory: {input_dir}")

    if not genome_path.exists():
        raise FileNotFoundError(f"Genome file not found: {genome_path}")
    if not genome_path.is_file():
        raise ValueError(f"Genome file path is not a file: {genome_path}")

    if main_output_dir.exists() and main_output_dir.is_file():
        raise ValueError(f"Output path exists as a file, not a directory: {main_output_dir}")

    bdg_metapeaks_dir = main_output_dir / "bdg_metapeaks"
    bdg_binary_matrices_dir = main_output_dir / "bdg_binary_matrices"
    bdg_bedgraphs_dir = main_output_dir / "bdg_bedGraphs"

    bdg_metapeaks_dir.mkdir(parents=True, exist_ok=True)
    bdg_binary_matrices_dir.mkdir(parents=True, exist_ok=True)
    bdg_bedgraphs_dir.mkdir(parents=True, exist_ok=True)

    base_name = input_dir.name

    concatenated_bed = main_output_dir / f"{base_name}_bedcat.bed"
    sorted_bed = main_output_dir / f"{base_name}_bedcat.sorted.bed"

    bedgraph_output_file = bdg_bedgraphs_dir / f"{base_name}_bedcat.bedGraph"
    macs2_output_file = bdg_metapeaks_dir / f"{base_name}_bedcat_metapeaks.bed"

    interval_files, skipped_files = find_bed_like_files(input_dir)

    if skipped_files:
        print("Skipped non-BED-like files:")
        for skipped_file, reason in skipped_files:
            print(f"  {skipped_file.name}: {reason}")

    if not interval_files:
        raise ValueError(f"No BED-like interval files found in input directory: {input_dir}")

    print(f"Found {len(interval_files)} BED-like files in: {input_dir}")
    for interval_file in interval_files:
        print(f"  accepted: {interval_file.name}")

    record_count = concatenate_interval_files(interval_files, concatenated_bed)
    if record_count == 0:
        raise ValueError(f"No interval records were found in the input directory: {input_dir}")

    print(f"Temporary concatenated BED written to: {concatenated_bed}")
    print(f"Concatenated BED record count: {record_count}")

    sort_cmd = [
        "sort",
        "-k1,1",
        "-k2,2n",
        str(concatenated_bed),
    ]

    sort_env = dict(os.environ)
    sort_env["LC_ALL"] = "C"

    print("Running command:", " ".join(sort_cmd))

    with sorted_bed.open("w") as out_handle:
        sort_result = subprocess.run(
            sort_cmd,
            stdout=out_handle,
            stderr=subprocess.PIPE,
            text=True,
            env=sort_env,
        )

    if sort_result.returncode != 0:
        raise RuntimeError(
            "UNIX sort failed.\n"
            f"Command: {' '.join(sort_cmd)}\n"
            f"stderr:\n{sort_result.stderr}"
        )

    if not sorted_bed.exists():
        raise RuntimeError(f"Sorted BED file was not created: {sorted_bed}")

    if sorted_bed.stat().st_size == 0:
        raise RuntimeError(
            "Sorted BED file was created but is empty.\n"
            f"Inspect these files:\n  {concatenated_bed}\n  {sorted_bed}"
        )

    print(f"Temporary sorted BED written to: {sorted_bed}")
    print(f"Sorted BED file size (bytes): {sorted_bed.stat().st_size}")

    genome_chroms = load_genome_chroms(genome_path)
    genome_chrom_set = set(genome_chroms)
    bed_chrom_set = load_bed_chroms(sorted_bed)
    shared_chroms = bed_chrom_set & genome_chrom_set

    print(f"Genome file chromosomes: {len(genome_chrom_set)}")
    print(f"Sorted BED chromosomes: {len(bed_chrom_set)}")
    print(f"Shared chromosomes: {len(shared_chroms)}")

    if bed_chrom_set:
        print("Example BED chromosomes:", ", ".join(sorted(list(bed_chrom_set))[:10]))
    if genome_chrom_set:
        print("Example genome-file chromosomes:", ", ".join(genome_chroms[:10]))

    if len(shared_chroms) == 0:
        raise RuntimeError(
            "No chromosome names are shared between the sorted BED file and the genome file.\n"
            "This usually means a naming mismatch such as 'chr1' versus '1'.\n"
            f"Inspect these files:\n  {sorted_bed}\n  {genome_path}"
        )

    genomecov_cmd = [
        "bedtools", "genomecov",
        "-bga",
        "-g", str(genome_path),
        "-i", str(sorted_bed),
    ]

    print("Running command:", " ".join(genomecov_cmd))

    with bedgraph_output_file.open("w") as out_handle:
        genomecov_result = subprocess.run(
            genomecov_cmd,
            stdout=out_handle,
            stderr=subprocess.PIPE,
            text=True,
        )

    if genomecov_result.returncode != 0:
        raise RuntimeError(
            "bedtools genomecov failed.\n"
            f"Command: {' '.join(genomecov_cmd)}\n"
            f"stderr:\n{genomecov_result.stderr}"
        )

    if not bedgraph_output_file.exists():
        raise RuntimeError(f"BedGraph output file was not created: {bedgraph_output_file}")

    if bedgraph_output_file.stat().st_size == 0:
        raise RuntimeError(
            "BedGraph output file was created but is empty.\n"
            f"Inspect these files:\n  {sorted_bed}\n  {genome_path}\n  {bedgraph_output_file}"
        )

    print(f"BedGraph written to: {bedgraph_output_file}")
    print(f"BedGraph file size (bytes): {bedgraph_output_file.stat().st_size}")

    cmd = [
        "macs2", "bdgpeakcall",
        "-i", str(bedgraph_output_file),
        "-o", str(macs2_output_file),
        "-c", str(args.cutoff),
        "-l", str(args.min_length),
        "-g", str(args.max_gap),
    ]

    print("Running command:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    if not macs2_output_file.exists():
        raise RuntimeError(f"MACS2 finished but output file was not created: {macs2_output_file}")

    print(f"MACS2 metapeaks written to: {macs2_output_file}")
    print("Temporary files kept for debugging:")
    print(f"  {concatenated_bed}")
    print(f"  {sorted_bed}")

    try:
        pybedtools.cleanup()
    except Exception:
        pass


if __name__ == "__main__":
    main()

