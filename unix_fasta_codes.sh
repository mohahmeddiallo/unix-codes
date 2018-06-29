#!/bin/bash
set -e
set -u
set -o pipefail


#Converts all bases to upper case letters
sed '/^[^>]/ y/acgt/ACGT/' input.fasta > output.fasta
#'/^[^>]/ This option excludes any line starting with '>'
#y transcribes the argument after first '/' into argument after second '/'



#This reverse transcribes rna to dna
sed '/^[^>]/ y/uU/tT/' uracil.fasta > thymine.fasta


#Removes ambiguous bases
sed '/^[^>]/ s/n//g' input.fasta > output.fasta


#Search through info line of fasta file for a regular pattern and
#produce a file with only those info with the pattern and their sequences-using 'grep' and 'sed'
grep -A1 "Nematoda" Database.fasta | sed '/^-/ s/--//'| sed '/^$/d' > Nematoda.fasta
#grep -A This command prints the line with the regular expression (info line) and the next one after it (seq line)
#first 'sed' Command to remove the third line with -- leaving an empty line
#second 'sed' command removes empty lines


#Converting from multiline fasta file to single line fasta file
perl -pe '/^>/ ? print "\n" : chomp' SILVA_128_LSURef_tax_silva_dna.fasta > SILVA_128_LSURef_tax_silva_dna_sl.fasta

#Remove first letter of a string
sed 's/^.\{1\}//g' gb203_pr2_all_10_28_97p_uparse.txt > gb203_pr2_all_10_28_97p_uparse_1.txt

