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

# The awk way
awk '/^>/ {printf("\n%s\n",$0)}
    /^>/ { printf("%s",$0);}
    END {printf("\n");}' < file.fa
# Another awk way
awk '/^>/ {printf("\n%s\n",$0);next; }
    { printf("%s",$0);}
    END {printf("\n");}' < file.fa
    

#Remove first letter of a string
sed 's/^.\{1\}//g' gb203_pr2_all_10_28_97p_uparse.txt > gb203_pr2_all_10_28_97p_uparse_1.txt

#Converting between alignment formats
seqmagick convert --output-format nexus --alphabet dna input.fasta output.nex


# For using file names as description lines for sequence entries in single entry .fasta files
for f in *.fasta; do name=$(sed 's/.fasta//' <<< $f); sed "/^>/ s/.*/>"$name"/" $f >> 18S_09_2019.fasta; done

for f in *.fasta; do name=$(sed 's/.fasta//' <<< $f); sed "/^>/ s/.*/>"$name"/" $f > ../"$name"_nrm.fasta; done

# For morphometrics
awk 'BEGIN{FS="\t"; print "Char", "\t", "Average","\t", "Range"} \
    {sum=0; for(i=2; i<=NF; i++) sum=sum+$i} \
    {max=0; for(i=2; i<=NF; i++) if(max<$i) max=$i} \
    {min=$2; for(i=2; i<=NF; i++) if(min>$i) min=$i} \
    {print $1, "\t", sum/(NF-1), "\t", min"–"max}' males_morphometrics.txt |column -t

# read lengths in fastq files for length distribution

awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' file.fastq


# Removes records with the PATTERN as id. It creates a backup of the input fastq file
sed -i.bak '/^@HWI-D00466:116:CC62WANXX:3:1102:7363:63646 1:N:0:GCACACG/,+3d' reads.fastq

# Remove short sequences [In this case shorter than 50]
grep -B 1 "[^>].\{50,\}" all.otus.txt


# Generate consensus sequences using awk
for iter in `seq $(awk 'BEGIN{FS=""} END{print NF}' $1)`
do
structure+=$(awk -v field=$iter 'BEGIN { FS = "" }
{if(count[$(field)]++ >= max) max = count[$(field)]}
END {for ( i in count ) if(max == count[i]) print i }' $1)
done

echo $structure > $2
structure=""

# getting list of filenames separated by comma for trinity
for fn in *2.fasta; do printf $(cat <<< $fn),; done > right


#Appending sequentially increasing number to specific lines denoted by a regular expression

awk 'BEGIN{ c = 10001} /^>/ {print $0 c++} \
    /^[^>]/ {print $0}' Onchocerca_volvulus_PRJNA37761.faa > Onchocerca_volvulus_PRJNA37761_modified.faa


# Prints for each entry sequence descriptions line and lengths of sequence
awk '/^>/ {print $0; next; } {print length($0)}' Aphelenchus_avenae_PRJNA236622_Trinity_first.fasta


# Takes a list of expressions line by line from a file
# Looks for them in a file with the patters and their
# associated string separated by a semicolon (eg. abbr and it its meaning)
# Takes both and performs substitution of abbreviation with the meanings wherever found in a different file.
while read line; do expr=$line; subs=$(grep $expr  NemaPhyloGenomics.csv | cut -d';' -f1 | cut -d' ' -f1,2) ; echo $subs; ; gsed -i "s/"$expr"/"$subs"/g" "spirurina_consensus.tre"; done < abrev
 

# This is the correct way of doing [ grep -v -A 1 '>Mel' inputfile ]
grep -n -A 1 '>Mel' OrthoGroup1438.fa |gsed -n 's/^\([0-9]\{1,\}\).*/\1d/p'| gsed -f - OrthoGroup1438.fa

 
# Trimming the ends of alignments using Gblocks

for FileName in *.fasta
do
sequences=`grep -c \> $FileName`
half_of_sequences=`expr $sequences / 2`
b1=`expr $half_of_sequences + 1` #Sets b1 as the floating point calculation (# of sequences)/2 + 1. If (# of sequences)/2 is not a whole number, the output is truncated so that is why 1 is added to each value.
b2=$b1 #Sets b2 = b1
b3=50 #default is 8
b4=10 #default is 10
b5=a #default is none
Gblocks $FileName -t=p -p=n -b1=$b1 -b2=$b2 -b3=$b3 -b4=$b4 -b5=$b5 # Executes Gblocks


rename 's/.fasta-gb/_cleaned.fasta/g' $FileName-gb
done


#Removing underscores '_' from .svg tree files

gsed -E '/^\s{6}>\w*/ s/_/ /g' input.svg > input_edit.svg
 
 
#Counting the number of bases (including multiline fasta formats)
awk '/^[^>]/ {printf("%s", $0)} END{printf("\n")}' abi/2763_1.fasta | wc -c

#Even simpler
awk '/^[^>]/' abi/2763_1.fasta | wc -c


#NCBI taxonomy to usearch
(^>.*)\.(\d) (.*) (.+) small.*
\1_\2;tax=g:\3,s:\4;


#PHASE OUTPUT
awk -F'[),]' '{for(i=1;i<=NF;i++){ if($i ~ /^:[0-9]{1,3}\./){print $i} } }' ENOPL_211005_exV1V9_nSS_NNNN_REV_RNA16A_MCMC.cnt > outputfile.txt



# MULTITHREADING
for x in *.fa
do
while [ $(ps -Af | grep "blast" | wc -l) -gt 12 ]
do
sleep 5
done

blast $x &...
sleep 1

done



# Converting Dieter's database to usearch format
awk -F'[ ;]' '/[^>]/ {accession=$1; domain=$2;phylum=$3;class=$4;order=$5;family=$6;genus=$7;species=$8; getline; printf("%s.%s_U;tax=d:%s,p:%s,c:%s,o:%s,f:%s,g:%s,s:%s_%s\n%s\n", accession,length,domain,phylum,class,order,family,genus,genus,species, $0)}' Nematode_ref_v10022022_final.fasta > Nematode_ref_v10022022_final_SINTAX.fasta
