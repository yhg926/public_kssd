<div align="center">
  <img src="Kssd_Logo" alt="KSSD logo" width="180">

# KSSD

**K-mer substring space sampling / shuffling decomposition for large-scale DNA sequence sketching, resemblance estimation, and containment analysis**

</div>

KSSD is a Linux command-line tool for sketching large collections of DNA sequences and estimating sequence resemblance or containment from reduced k-mer representations. It accepts FASTA and FASTQ input, including gzipped files, and can also sketch streamed sequence data, for example from SRA Toolkit commands.

KSSD differs from ordinary per-dataset k-mer sampling methods by using a predefined shuffled k-mer substring space. References and queries are mapped into the same sampled k-mer subspace, enabling efficient large-scale comparison while preserving both resemblance- and containment-oriented signals.

If you use KSSD, please cite the Genome Biology paper listed in [Citation](#citation).

---

## Contents

- [Features](#features)
- [Installation](#installation)
- [Command overview](#command-overview)
- [Quick tutorial](#quick-tutorial)
- [Core workflow](#core-workflow)
- [Choosing k-mer and reduction parameters](#choosing-k-mer-and-reduction-parameters)
- [Generate a shuffled k-mer substring space](#generate-a-shuffled-k-mer-substring-space)
- [Sketch references](#sketch-references)
- [Sketch queries](#sketch-queries)
- [Compare sketches](#compare-sketches)
- [Sketch from streaming input](#sketch-from-streaming-input)
- [Combine query batches](#combine-query-batches)
- [Set operations](#set-operations)
- [Reverse sketches to k-mer sets](#reverse-sketches-to-k-mer-sets)
- [Output: distance.out](#output-distanceout)
- [Notes and limitations](#notes-and-limitations)
- [Citation](#citation)
- [License](#license)

---

## Features

- Sketch DNA sequences from FASTA or FASTQ files.
- Supports gzipped and uncompressed input files.
- Builds indexed reference sketch databases.
- Searches query sketches against reference sketch databases.
- Computes pairwise resemblance and containment statistics.
- Supports Jaccard-like and containment-oriented output modes.
- Supports sketch union, intersection, subtraction, unique union, and pan-sketch combination operations.
- Supports streaming input through `--pipecmd`.
- Supports multi-threading when compiled with OpenMP.
- Designed for Linux systems.

---

## Installation

### Requirements

KSSD is currently intended for Linux. A typical build requires:

- `gcc`
- `make`
- `zlib`
- OpenMP support, usually available through GCC

On Ubuntu or Debian:

```bash
sudo apt-get update
sudo apt-get install -y build-essential zlib1g-dev
```

### Build from source

```bash
git clone https://github.com/yhg926/public_kssd.git
cd public_kssd
make
```

After compilation, the executable should be available as:

```bash
./kssd
```

You can optionally add the repository directory to your `PATH`.

---

## Command overview

The executable exposes the following top-level subcommands:

| Command | Purpose |
|---|---|
| `kssd shuffle` | Generate a shuffled k-mer substring-space file (`.shuf`). |
| `kssd dist` | Sketch sequences, build/index reference sketches, and estimate distances or containment. |
| `kssd set` | Perform set operations on sketches. |
| `kssd reverse` | Recover k-mer sets from KSSD sketches. |
| `kssd composite` | Advanced metagenomic composition analysis command. |

Most users will mainly use `shuffle`, `dist`, and `set`.

---

## Quick tutorial

The repository includes a small example dataset under `test_fna/`.

```bash
cd test_fna

# 1. Sketch and index reference sequences.
../kssd dist -L ../shuf_file/L3K10.shuf -o reference ./seqs1
../kssd dist -o reference reference

# 2. Sketch query sequences with the same .shuf file.
../kssd dist -L ../shuf_file/L3K10.shuf -o query ./seqs2

# 3. Search queries against the reference database.
../kssd dist -r reference -o distout query

# 4. Or compute pairwise distances among the reference sketches.
../kssd dist -r reference -o distout reference
```

The main output file is:

```text
distout/distance.out
```

---

## Core workflow

A common KSSD analysis has four stages:

```text
Input FASTA/FASTQ files
        |
        v
Sketch sequences with a .shuf file
        |
        v
Build/index reference sketches
        |
        v
Compare query sketches against references
        |
        v
distance.out
```

Important: references and queries must be sketched using the same `.shuf` file. Otherwise, their sampled k-mer spaces are not comparable.

---

## Choosing k-mer and reduction parameters

KSSD uses half-length notation for k-mers.

For example:

```bash
-k 8
```

means a full k-mer length of `16`.

### Common starting points

| Data type | Suggested `-k` | Full k-mer length | Suggested reduction level |
|---|---:|---:|---:|
| Bacterial / prokaryotic genomes | `8` | 16 | `-L 2` or `-L 3` |
| Intermediate-size genomes | `9` | 18 | `-L 3` or `-L 4` |
| Mammalian genomes or large metagenomic datasets | `10` or `11` | 20 or 22 | `-L 4` or `-L 5` |

The expected dimensionality reduction rate for level `L` is approximately:

```text
16^L
```

Examples:

| `L` | Approximate reduction |
|---:|---:|
| 2 | 256× |
| 3 | 4,096× |
| 4 | 65,536× |
| 5 | 1,048,576× |

Higher reduction levels produce smaller sketches and faster comparisons, but may reduce sensitivity.

---

## Generate a shuffled k-mer substring space

A `.shuf` file defines the shuffled k-mer substring space used for sketching.

```bash
kssd shuffle \
  -k <half_kmer_length> \
  -s <half_substring_length> \
  -l <reduction_level> \
  -o <output_prefix_or_file>
```

Example:

```bash
kssd shuffle -k 8 -s 5 -l 3 -o bacteria_L3K8.shuf
```

Options:

| Option | Meaning | Default in command help |
|---|---|---:|
| `-k` | Half-length of the k-mer. Full k-mer length is `2 × k`. | `8` |
| `-s` | Half-length of the k-mer substring. | `5` |
| `-l` | Dimensionality reduction level. Expected reduction is `16^l`. | `2` |
| `-o` | Output prefix or file name for the generated `.shuf` file. | `./default` |
| `--usedefault` | Use the default prokaryote-oriented setting: `k=8`, `s=5`, `l=2`. | off |

Parameter relationship:

```text
l < s < k
```

You can also use the provided `.shuf` files under `shuf_file/`.

---

## Sketch references

There are two common ways to sketch references.

### A. Use an existing `.shuf` file

```bash
kssd dist \
  -L <shuf_file> \
  -o <reference_output_dir> \
  <reference_sequence_dir>
```

Example:

```bash
kssd dist -L shuf_file/L3K10.shuf -o reference genomes/
```

Then build/index the reference sketch database:

```bash
kssd dist -o reference reference
```

### B. Let `kssd dist` generate a default `.shuf` file

```bash
kssd dist \
  -L <reduction_level> \
  -k <half_kmer_length> \
  -o <reference_output_dir> \
  <reference_sequence_dir>
```

Example:

```bash
kssd dist -L 3 -k 8 -o reference genomes/
```

When `-L` receives an integer instead of an existing `.shuf` file, KSSD generates and uses a shuffled space internally in the output directory.

---

## Sketch queries

Queries must be sketched with the same `.shuf` file used for the references.

```bash
kssd dist \
  -L <same_shuf_file_used_for_references> \
  -o <query_output_dir> \
  <query_sequence_dir>
```

Example:

```bash
kssd dist -L shuf_file/L3K10.shuf -o query query_fastq/
```

---

## Compare sketches

### Search queries against references

```bash
kssd dist \
  -r <reference_sketch_dir> \
  -o <output_dir> \
  <query_sketch_dir>
```

Example:

```bash
kssd dist -r reference -o search_out query
```

The main result is:

```text
search_out/distance.out
```

### Pairwise comparison among references

```bash
kssd dist \
  -r <reference_sketch_dir> \
  -o <output_dir> \
  <reference_sketch_dir>
```

Example:

```bash
kssd dist -r reference -o pairwise_out reference
```

### Useful `dist` options

| Option | Meaning | Default / notes |
|---|---|---|
| `-k, --halfKmerlength INT` | Half k-mer length. Full k-mer length is `2 × k`. | `8` |
| `-p, --threadN INT` | Number of threads. | all available threads if OpenMP is enabled |
| `-l, --list file` | File containing paths to query sequences. | optional |
| `-L, --DimRdcLevel INT/file` | Reduction level or `.shuf` file. | `2` |
| `-m, --maxMemory NUM` | Maximum memory in GB. | system memory |
| `-n, --LstKmerOcrs INT` | Minimum k-mer occurrence threshold in FASTQ input. | `1`, capped internally at `7` |
| `-Q, --quality INT` | Filter k-mers with lowest base quality below this Phred value. | `0` |
| `-r, --reference_dir path` | Reference sketch/database directory. | optional |
| `-o, --outdir path` | Output directory. | `.` |
| `-N, --neighborN_max INT` | Maximum number of nearest reference genomes to report. | command help says `1`; internal default is `0` |
| `-D, --mutDist_max FLT` | Maximum mutation distance allowed for output. | `1` |
| `-M, --metric 0/1` | Output metric: `0` = Jaccard, `1` = containment. | `0` |
| `-O, --outfields 0/1/2` | Output fields: distance, q-values, confidence intervals. Later levels include earlier ones. | `2` |
| `--correction 0/1` | Correct shared k-mer counts or not. | `0` |
| `-A, --abundance` | Abundance-estimation mode. | off |
| `-u, --dedup` | Ignore repeated k-mers in reference. | off |
| `--keepcofile` | Keep intermediate `.co` files. | off |
| `--pipecmd cmd` | Read input from a pipe command. | optional |
| `--keepskf` | Keep shared-kmer-count file. | off |
| `-f, --skf path` | Shared-kmer-count file path. | optional |
| `--byread` | Sketch input by read. | off |

---

## Sketch from streaming input

KSSD can sketch sequence data from a command pipeline. This is useful for SRA accessions.

Example using `fastq-dump`:

```bash
kssd dist \
  -L <shuf_file> \
  -n 2 \
  -o <output_dir> \
  --pipecmd "fastq-dump --skip-technical --split-spot -Z" \
  ERR000001
```

Or after prefetching the SRA file:

```bash
prefetch ERR000001

kssd dist \
  -L <shuf_file> \
  -n 2 \
  -o <output_dir> \
  --pipecmd "fastq-dump --skip-technical --split-spot -Z" \
  /path/to/ERR000001.sra
```

---

## Combine query batches

If query sketches were generated in multiple batches using the same `.shuf` file, they can be combined:

```bash
kssd dist \
  -o <combined_output_dir> \
  <query_batch_1> <query_batch_2> ...
```

Example:

```bash
kssd dist -o combined_queries query_batch1 query_batch2 query_batch3
```

All batches must have been generated with compatible sketch parameters.

---

## Set operations

`kssd set` operates on combined sketch directories.

### Union

```bash
kssd set -u -o <union_output_dir> <combined_sketch_dir>
```

Example:

```bash
kssd set -u -o union_out query/qry
```

### Unique union

```bash
kssd set -q -o <unique_union_output_dir> <combined_sketch_dir>
```

### Intersection with a pan-sketch

```bash
kssd set -i <pan_sketch_dir> -o <intersection_output_dir> <combined_sketch_dir>
```

Example:

```bash
kssd set -i union_out -o intersect_out query/qry
```

### Subtraction of a pan-sketch

```bash
kssd set -s <pan_sketch_dir> -o <subtraction_output_dir> <combined_sketch_dir>
```

Example:

```bash
kssd set -s union_out -o subtract_out query/qry
```

### Other `set` options

| Option | Meaning |
|---|---|
| `-c, --combin_pan` | Combine pan files into a combined sketch file. |
| `-p, --threads INT` | Number of threads. |
| `-P, --print` | Print genome names. |
| `-g, --grouping file.tsv` | Group genomes by an input category file. |
| `-o, --outdir path` | Output directory. |

---

## Reverse sketches to k-mer sets

`kssd reverse` is an advanced command for recovering k-mer sets from KSSD sketch directories.

```bash
kssd reverse \
  -L <shuf_file> \
  -o <output_dir> \
  <co_dir>
```

Options:

| Option | Meaning |
|---|---|
| `-L, --shufFile path` | Provide the `.shuf` file used for sketching. |
| `-o, --outdir path` | Output directory for recovered k-mer files. |
| `-p, --threads INT` | Number of threads. |
| `-b, --byreads` | Recover k-mers from sketches generated by read. |

---

## Output: distance.out

The main comparison result is written to:

```text
<output_dir>/distance.out
```

Typical columns include:

| Column | Description |
|---|---|
| `Qry` | Query sequence or sketch name. |
| `Ref` | Reference sequence or sketch name. |
| `Shared_k` | Number of shared sampled k-mers between query and reference sketches. |
| `Ref_s` | Reference sketch size. |
| `Qry_s` | Query sketch size. |
| `Jaccard` | Jaccard coefficient. |
| `MashD` | Mash-style distance. |
| `ContainmentM` | Containment measurement. |
| `AafD` | AAF-style distance. |
| `Jaccard_CI` | 95% confidence interval for Jaccard. |
| `MashD_CI` | 95% confidence interval for Mash distance. |
| `ContainmentM_CI` | 95% confidence interval for containment measurement. |
| `AafD_CI` | 95% confidence interval for AAF distance. |
| `P-value(J)` | P-value for Jaccard. |
| `P-value(C)` | P-value for containment. |
| `FDR(J)` | False discovery rate for Jaccard. |
| `FDR(C)` | False discovery rate for containment. |

The exact fields may depend on options such as `-M` and `-O`.

---

## Notes and limitations

- KSSD currently targets Linux. macOS and Windows are not officially supported.
- References and queries must be sketched with the same `.shuf` file.
- Report the `.shuf` file, `-k`, `-s`, reduction level, KSSD version, and command lines in reproducible analyses.
- Very high reduction levels reduce sketch size and runtime but may reduce sensitivity.
- Some advanced commands, especially `composite`, are exposed by the executable but are less documented than the main `shuffle`, `dist`, and `set` workflow.

---

## Citation

If you use KSSD, please cite:

```text
Yi, H., Lin, Y., Lin, C. et al.
Kssd: sequence dimensionality reduction by k-mer substring space sampling enables real-time large-scale datasets analysis.
Genome Biology 22, 84 (2021).
https://doi.org/10.1186/s13059-021-02303-4
```

---

## License

KSSD is distributed under the Apache License, Version 2.0. See `LICENSE.txt` for details.

---

## Contact

Please use the GitHub Issues page for bug reports, questions, and feature requests:

```text
https://github.com/yhg926/public_kssd/issues
```
