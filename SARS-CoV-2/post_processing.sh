#!/usr/bin/env bash
DATE=$1

cd $DATE && find . -type f -name "*consensus.fa" -mindepth 2 -exec mv -i -- {} . \;
mv /home/vitalite/iSeq_share/$DATE/*.fastq.gz .
mkdir NML
mv *.consensus.fa NML
mv *.fastq.gz NML

echo "generating sample table and preparing files for transfer..."

tail -n +15 /home/vitalite/iSeq_share/$DATE/SampleSheet.csv | cut -f2 -d ',' | awk 'NR<2{print $0;next}{print $0| "sort -r"}' > tmp.txt
awk -vDATE="$DATE" 'BEGIN{FS=OFS="\t"}{print $0,( NR==1 ? "run_name" OFS "barcode" OFS "primer_scheme" OFS "analytical_protocol" : DATE OFS "NA" OFS "freed" OFS "SIGNAL_v1.4.3")}' tmp.txt > tmp2.txt
printf "consensus_filename\n$(ls NML/*consensus.fa)" | awk 'NR<2{print $0;next}{print $0| "sort -r"}' > tmp3.txt
paste -d '\t' tmp2.txt tmp3.txt | awk 'BEGIN{FS=OFS="\t"}{print $0,( NR==1 ? "GISAID_virus_name" OFS "GISAID_accession_ID" : "NA" OFS "NA")}' > NML_Moncton_$DATE.0.tsv
rm tmp.txt
rm tmp3.txt
rm tmp2.txt

mv NML_Moncton_$DATE.0.tsv NML
tail -q -n +2 lineage.assignments.tsv | cut -d '_' -f 2- > new.tsv
head -n 1 lineage.assignments.tsv > new2.tsv
cat new2.tsv new.tsv > lineage_assignments_$DATE.tsv
rm new.tsv
rm new2.tsv
rm lineage.assignments.tsv
mv summary.html summary_$DATE.html
mv freebayes.lineage.assignments.tsv freebayes.lineage.assignments_$DATE.tsv
mv /home/vitalite/covid-19-signal/archive/$DATE/plots/ncov_tools_graphs_depth_by_position.pdf /home/vitalite/covid-19-signal/archive/$DATE/plots/ncov_depth_position_$DATE.pdf
mv /home/vitalite/covid-19-signal/archive/$DATE/plots/ncov_tools_graphs_amplicon_covered_fraction.pdf /home/vitalite/covid-19-signal/archive/$DATE/plots/ncov_amplicon_covered_$DATE.pdf

#python /home/vitalite/NML_CanCOGeN_new-brunswick.py -u NML


echo "finished"






