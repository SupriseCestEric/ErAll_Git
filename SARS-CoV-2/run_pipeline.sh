#!/usr/bin/env bash

OUTDIR=$1

conda activate signal
snakemake -kp --configfile config.yaml --cores=2 --use-conda --conda-prefix=$PWD/.snakemake/conda all
sleep 30
snakemake -p --configfile config.yaml --cores=2 --use-conda --conda-prefix=$PWD/.snakemake/conda postprocess
sleep 30
cd results_dir && find . -name "*_*" -type f -print0 | while read -d $'\0' f; do mv "$f" "$(echo $f | sed 's/_/\./g')"; done &&
cd /home/vitalite/covid-19-signal
conda deactivate
sleep 10 &&
conda activate ncov-qc && cd /home/vitalite/covid-19-signal/ncov-tools &&
snakemake --cores=2 -s workflow/Snakefile all_qc_sequencing --configfile config.yaml &&
cd /home/vitalite/covid-19-signal &&
sleep 10 &&
mkdir archive/$OUTDIR && mv -v results_dir/* archive/$OUTDIR &&
mv /home/vitalite/covid-19-signal/ncov-tools/plots /home/vitalite/covid-19-signal/archive/$OUTDIR &&
rm -rf /home/vitalite/covid-19-signal/ncov-tools/qc_sequencing &&
conda deactivate
