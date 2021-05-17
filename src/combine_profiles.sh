#!/usr/bin/bash

parent_dir="/project/shefflab/processed/pepatac/meta_cancer/"
work_dir="/project/shefflab/processed/pepatac/meta_cancer/profiling/"

sample_list="${work_dir}sample_list.txt"

function build_profile {
  profile_file="${parent_dir}/results_pipeline/${1}/PEPATAC_profile.tsv"
  cp "${profile_file}" "${workdir}profile.tmp"
  grep 'File_mb' "${parent_dir}/results_pipeline/${1}/PEPATAC_log.md" | awk '{print $3}' > "${workdir}file_mb.tmp"
  grep 'Raw_reads' "${parent_dir}/results_pipeline/${1}/PEPATAC_log.md" | awk '{print $3}' > "${workdir}raw_reads.tmp"
  echo $1 > "${workdir}sample.tmp"
  awk 'FNR > 2 {print $0}' "${workdir}profile.tmp" | head -n 1 | awk -v OFS='\t' '{print $0,"runtime_seconds","file_mb","raw_reads","sample"}' > "${workdir}header.tmp"
  sed -i '1,3d' "${workdir}profile.tmp"
  awk '{print $4}' "${workdir}profile.tmp" | awk -F: '{ print ($1 * 3600) + ($2 * 60) + $3 }' > "${workdir}seconds.tmp"
  paste "${workdir}profile.tmp" "${workdir}seconds.tmp" "${workdir}file_mb.tmp" "${workdir}raw_reads.tmp" "${workdir}sample.tmp" > "${workdir}profile_${1}.tmp"
  cat "${workdir}header.tmp" "${workdir}profile_${1}.tmp" > "${1}_profile.tsv"
  find $workdir -type f -name '*.tmp' -exec rm {} \;
}

while read sample; do
  echo $sample
  build_profile $sample
done <$sample_list