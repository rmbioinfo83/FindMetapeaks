# FindMetapeaks

FindMetapeaks builds metapeaks from a directory of BED-like peak files, then generates a **sample-by-metapeak count matrix**.

The current pipeline does the following:

1. Detects BED-like interval files in the input directory by **content**, not filename extension.
2. Concatenates those files into one combined BED file.
3. Sorts the combined BED file.
4. Uses `bedtools genomecov -bga` to generate a BedGraph coverage track.
5. Uses `macs2 bdgpeakcall` to call metapeaks from that BedGraph.
6. Uses `bedtools intersect -c` to count, for each sample file, how many peaks overlap each metapeak.
7. Writes a **count matrix** with:
   - **rows = samples**
   - **columns = metapeaks**
   - first column name = `sample_filenames`

## Requirements

Validated requirements:

- **Python 3.10 or 3.11 recommended**
- `macs2` installed and available on `PATH`
- `bedtools` installed and available on `PATH`
- standard Unix `sort` available on `PATH`
- `pybedtools` installed in the same Python environment

Notes:

- In testing, **Python 3.12 was not a good fit for MACS2**.
- A clean install was successfully validated in **WSL Ubuntu 22.04** with Python 3.11.
- On **WSL**, install and build FindMetapeaks from the **Linux filesystem** (for example `~/FindMetapeaks`), **not** from a mounted Windows path such as `/mnt/c/...`.

## Installation

### 1. Create and activate a Python environment

Example with Python 3.11:

```bash
python3.11 -m venv findmetapeaks_test_env
source findmetapeaks_test_env/bin/activate
python -m pip install --upgrade pip
```

### 2. Install Python dependencies

```bash
python -m pip install macs2 pybedtools
```

### 3. Make sure external tools are available

`bedtools` must already be installed and on `PATH`.

Quick checks:

```bash
which macs2
macs2 --version
which bedtools
bedtools --version
python -c "import pybedtools; print(pybedtools.__version__)"
```

### 4. Install FindMetapeaks from a local checkout

From the root of the package directory (the directory containing `pyproject.toml`):

```bash
python -m pip install .
```

For development installs, you can also use:

```bash
python -m pip install -e .
```

### 5. Install FindMetapeaks from GitHub

After the repository is published on GitHub, users can install directly from the repository with:

```bash
python -m pip install git+https://github.com/<YOUR_GITHUB_USERNAME>/FindMetapeaks.git
```

If you later tag releases, you can also install a specific tag or branch, for example:

```bash
python -m pip install git+https://github.com/<YOUR_GITHUB_USERNAME>/FindMetapeaks.git@v0.1.0
```

### 6. Confirm the command was installed

```bash
bdgFindMetapeaks --help
```

## Input requirements

### Input directory (`-i` / `--input`)

The input must be a directory containing **BED-like interval files**.

Important:

- File extensions do **not** have to be `.bed`
- `.narrowPeak` and other BED-like files are accepted
- Files are accepted by checking that they contain valid genomic interval records
- BED-like records must have at least the first 3 columns:
  - chromosome
  - start
  - end

Example:

```text
/project/.../bedcat_sub10000_CUTOFF/H3K27ac
```

### Genome file (`-G` / `--genome-file`)

This must be a standard 2-column genome file matching the genome build used for the input interval files.

Example format:

```text
chr1    248956422
chr2    242193529
chr3    198295559
```

### Output directory (`-o` / `--output-dir`)

This is the **main output directory** created by FindMetapeaks.

## Usage

```bash
bdgFindMetapeaks \
    -i /path/to/input_peak_directory \
    -o /path/to/main_output_dir \
    -G /path/to/genome_file \
    --cutoff 5.0 \
    --min-length 200 \
    --max-gap 30
```

Default parameters:

- `--cutoff 5.0`
- `--min-length 200`
- `--max-gap 30`

These defaults were chosen to match the MACS2 `bdgpeakcall` defaults used during development.

## Output structure

Given an input directory named `H3K27ac`, FindMetapeaks creates:

```text
main_output_dir/
├── H3K27ac_bedcat.bed
├── H3K27ac_bedcat.sorted.bed
├── bdg_bedGraphs/
│   └── H3K27ac_bedcat.bedGraph
├── bdg_metapeaks/
│   └── H3K27ac_bedcat_metapeaks.bed
└── bdg_count_matrices/
    └── H3K27ac_bedcat_metapeak_count_matrix.tsv
```

### Notes on outputs

- `H3K27ac_bedcat.bed` is the concatenated input BED file
- `H3K27ac_bedcat.sorted.bed` is the sorted version of that BED file
- these two files are currently **kept for debugging**
- the BedGraph file is kept in `bdg_bedGraphs/`
- the called metapeaks are kept in `bdg_metapeaks/`
- the non-binary count matrix is kept in `bdg_count_matrices/`

## Count matrix format

The saved matrix is **not binarized**.

It is a **count matrix** where each value is the number of sample peaks overlapping a given metapeak.

### Matrix orientation

- **rows = samples**
- **columns = metapeaks**

### First column

The first column is explicitly named:

```text
sample_filenames
```

This column contains the original input filenames, including file extensions.

### Metapeak column names

Metapeak columns are named from genomic coordinates in the format:

```text
chrom_start_end
```

For example:

```text
chr1    10001   10542
```

becomes:

```text
chr1_10001_10542
```

### Why counts instead of binary values?

The matrix stores **absolute overlap counts** so users can choose their own downstream transformation, including binarization if desired.

## Method summary

At a high level, the pipeline is:

1. Identify BED-like files in the input directory
2. Concatenate them into one combined BED
3. Sort the combined BED
4. Generate a BedGraph coverage track with `bedtools genomecov -bga`
5. Call metapeaks with `macs2 bdgpeakcall`
6. For each sample file, compute overlap counts against the metapeaks using `bedtools intersect -c`
7. Assemble the final sample-by-metapeak count matrix

## Validated behavior

The current implementation has been validated in two ways:

1. The metapeak BED output matched direct MACS2-generated results.
2. The final count matrix matched a previously generated reference matrix exactly, including per-sample metapeak count vectors.

## Troubleshooting

### `macs2` fails to install

If you are using Python 3.12, try Python 3.10 or 3.11 instead.

### `bedtools` not found

Make sure `bedtools` is installed and on `PATH`:

```bash
which bedtools
bedtools --version
```

### `bdgFindMetapeaks` command not found after install

Make sure:

- the correct virtual environment is activated
- the package install completed successfully

Then check:

```bash
which bdgFindMetapeaks
```

### WSL build or install fails under `/mnt/c/...`

If you are using WSL, copy the repository into the Linux filesystem (for example `~/FindMetapeaks`) and install it there instead of building from a mounted Windows path.
