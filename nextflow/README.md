# Nextflow scripts for processing the manuscript data

**NOTE**: You will need access to the BAM files (EGA study EGAS00001004670) to rerun the analysis.

## Prerequisites

General requirements:

- [nextflow](https://nextflow.io)
- A working installation of ichorCNA
- `readCounter` from https://github.com/shahcompbio/hmmcopy_utils

For `bam_size_select.py`:

- [arghandler](https://github.com/druths/arghandler)
- [joblib](https://joblib.readthedocs.io/en/latest/)
- [pysam](https://pysam.readthedocs.io/en/latest/index.html)
- [natsort](https://pypi.org/project/natsort/)

## Setup

Prepare the `samples.csv` file according to the template, or as follows:

- First column (`patient`): patient
- Second column (`samplename`): sample name
- Third column (`bamfile`): full path to the sample's BAM file
- Fourth column (`type`): either `tissue` or `cfdna` depending on the type

Copy the `bam_size_select.py` script somewhere in your PATH or where the analysis can find it
and put its path into `bam_size_select` in `nextflow.config`.

Adjust `gcwig`, `mapwig`, `centromere` and `ichorcna_script` to your ichorCNA installation (add also `ichorcna_libdir` if you're using a git checkout), and supply the full path to the panel of normals in `normal_panel`. Adjust the `profiles` section according to your setup (cluster, Docker, etc.).

## Running

Run it via nextflow:

```
nextflow run -c nextflow.config -profile <your profile> ichorcna.nf --index <path to samples.csv> --outdir <destination folder>
```

## Results

The script will perform BAM size selection if plasma samples are detected (`type` is `ctdna` in the CSV file) and will create individual folders for each patient. A summarized CSV file will be created in the root of the output directory, with summarized tumor fraction estimates.

