#!/bin/bash

for assembly in $(ls *.fasta)
do
	for CHR in $(cat /home/dennu/Desktop/ASSEMBLY/Chr_name.tsv | cut -f 1)
	do
		chr_new=$(grep $CHR /home/dennu/Desktop/ASSEMBLY/Chr_name.tsv | cut -f 2)
		sed -i -E "s/.*$CHR*\t\S*/$chr_new\t/g" $assembly
	done
done
