#!/usr/bin/env bash

DATE=$1
runName=$2
runID="$(bs -c canada list runs | grep $runName | cut -f4 -d' ')"

echo "Finding biosamples and raw files from the provided run ID ..."
mkdir $DATE
bs -c canada download run -i $runID --extension=csv -o $DATE
tail -n +16 $DATE/SampleSheet.csv | cut -f1 -d ',' > searchterms.txt
bs -c canada list biosamples | grep -F -f searchterms.txt | awk '{print $4}' > biosamples.txt

echo "Downloading fastq files associated to given samples ..."
while read x; do
  bs -c canada download biosample -i $x --extension=fastq.gz -o $DATE
done < biosamples.txt
rm searchterms.txt
mv biosamples.txt $DATE

cd $DATE
find . -type f -mindepth 2 -exec mv -i -- {} . \;
