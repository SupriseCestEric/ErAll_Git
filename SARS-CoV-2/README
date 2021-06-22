# Using the SIGNAL pipeline on Vitalite server

## Fetching files from basespace

Navigate to the </home/vitalite/iSeq_share> directory and execute the following command

bs -c canada list runs

The name of the run should resemble nCov-month-day or something similar. 
run the get_bs_files script with two arguments, the date and run name

./get_bs_files.sh YYYYMMDD nCov-month-day

This automatically downloads all the biosamples and fastq files associated with the run.
Afterwards, navigate to /home/vitalite/covid-19-signal and generate the sample table using

./generate_sample_table.sh -d /home/vitalite/iSeq_share/YYYYMMDD

This should automatically return "success!" and auto-fills the .csv used to point SIGNAL to the appropriate files. 

## Running SIGNAL, the post-processing script and ncov-tools

The wrapper script run_pipeline.sh handles all steps of the signal pipeline and plots with ncov-tools.
It needs to be run in iterperative moode (-i) because it will activate and use many conda environments in its own shell.

bash -i run-pipeline.sh YYYYMMDD

This stores all results in the results_dir folder, stores all plots in /home/vitalite/covid-19-signal/ncov-tools. 
However, it also automatically creates a YYYYMMDD directory in covid-19-signal/archive and moves all results there. 
Another script is used to timestamp results, find and transfer fastq / fasta files and auto-generate a TSV in NML requested format, and populates it. 
Be wary of the trailing backslash (YYMMDD/), which will not work if present. 

./post_processing.sh YYYYMMDD

## Transfering to NML

Finally, data is transfered to NML using the azure script provided by Winnipeg when all data has been revised for errors.

python /home/vitalite/NML_CanCOGeN_new-brunswick.py -u NML
