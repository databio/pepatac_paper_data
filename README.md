# PEPATAC paper analyses

## Paper figures

To produce the figures in the paper, you'll first need to download and process the samples found in the "paper_sample_table.csv" included file.

After downloading, you can process using the pipeline:
```
looper run paper_config.yaml
looper run paper_none_config.yaml
```

After completing, generate summary statistics:
```
looper report paper_config.yaml
looper report paper_none_config.yaml
```

This will produce output variants with prealignments and without for downstream comparisons.

The [included markdown file](src/PEPATAC_paper_plots.md) may be followed to reproduce the plots in R from the paper.

## Prealignment comparisons
To produce the prealignment timing comparison plots requires three primary steps.

### 1. Obtain source files

The mitochondrial (mtDNA) and human nuclear genome (hg38) aligning reads are originally derived from the following GEO accessions:
 - GSM2471255
 - GSM2471300
 - GSM2471249
 - GSM2471269
 - GSM2471245

### 2. Run source files through PEPATAC and *keep* prealignment BAM files

We want to extract mitochondrial reads, so we will keep all prealignment files. The default is to remove them to save disk space.  The included "source_library_config.yaml" is our PEP for these samples.

```
looper run source_library_config.yaml -x " --keep"
```

### 3. Extract mitochondrial and nuclear genome aligning reads

After these samples finish, we want to generate all of the various total read counts necessary of both mtDNA and hg38 aligning reads that we can combine in various ratios to generate 10-100% mixtures from 10M to 200M total reads per mixture.

`./generate_libraries.sh "mtDNA_reads" "hg38_reads" "/path/to/source_library_output/results_pipeline/"`

This is best accomplished on a cluster or a machine with upwards of 100GB of available RAM.

### 4. Analyze using prealignments and without

Set a environment variable that points to the directory containing your generated libraries named DATA.

`export DATA=/your/path/to/mtDNA_reads/`

Run each version of the PEP project using the same compute resources.

```
looper run prealignment_config.yaml --compute cores=8 mem=16000
looper run prealignment_none_config.yaml --compute cores=8 mem=16000
```

### 5. Produce comparison plots

```
Rscript PEPATAC_profile_aggregator.R /path/to/prealignment_config.yaml /path/to/prealignment_none_config.yaml $PROCESSED/pepatac/prealignment_comparison/yes/results_pipeline/ $PROCESSED/pepatac/prealignment_comparison/no/results_pipeline/ /path/to/your/output_dir
```