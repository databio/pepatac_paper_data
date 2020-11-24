# PEPATAC paper analyses

## Get the data files

The data used for the paper is available from public resources. Included in the `metadata/` subfolder is the 
"paper_sra_accessions.txt" file containing a list of sequence read archive accession numbers.

To obtain the files en masse, you can provide the entire file to NCBI's sra-tools' `fasterq-dump` function like so:
```
cat paper_sra_accessions.txt | xargs -n1 fasterq-dump -p -O /path/to/output_dir
```

To simplify use of downstream configuration files, you can also create an environment variable (`SRAFQ`) that points to this output directory containing your fastq files.

```
export SRAFQ=/path/to/output_dir
```

## Run the pipeline

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

The [included R markdown file](src/PEPATAC_paper_plots.Rmd) may be followed to reproduce the plots in R from the paper.

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