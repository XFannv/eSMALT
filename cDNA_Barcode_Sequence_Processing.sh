# generate bystrand ccs reads--------------------------------------------------------------

ccs -j 24 --by-strand --min-passes 3 --min-rq 0.95 --min-length 800 --max-length 1800 \
subreads.bam bystrand_ccs.bam 

samtools fastq bystrand_ccs.bam  >> ccs.fq

## mapping--------------------------------------------------------------

cat sample_list.txt | while read sample; do

minimap2 -cx map-pb --cs=long -t 12 -A 2 -B 2 -O 2,12 ref.mmi ${sample}.fq |
awk '{if($12==60 && $10>=900 && $8==0 && $9==1030) print $0}' |
cut -f1,24 > mapping/${sample}.cstag

done


### UMI Cluster~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# extract quality raw sequence
mkdir raw_qual_reads

cat cell_list.txt | cut -f1 | while read sample; do

cat mapping/${sample}.cstag | awk '{print"@"$1}' > raw_qual_reads/${sample}.id

cat demultiplex/${sample}.fq | paste - - - - | sed 's/ rc//g' |
awk 'NR==FNR{a[$1]=$1}NR>FNR{if($1 in a)print$0}' raw_qual_reads/${sample}.id - |
tr "\t" "\n" > raw_qual_reads/${sample}.fq

done


# extract quality clean sequence
mkdir clean_qual_reads

cat cell_list.txt | cut -f1 | while read sample; do

cat mapping/${sample}.cstag | cut -f2 |
sed -e 's/cs:Z://g' -e 's/=//g' |
sed -e 's/+[atcg]*//g' -e 's/-//g' |
sed -e 's/*ga/a/g' -e 's/*ct/t/g' |
sed -e 's/*a[tcg]/a/g' -e 's/*t[acg]/t/g' -e 's/*c[ag]/c/g' -e 's/*g[tc]/g/g' |
tr atcg ATCG | paste mapping/${sample}.cstag - | cut -f1,3 |
sed -e 's/^m/>m/g' | sed -e 's/\t/\n/g' > clean_qual_reads/${sample}.fa

done


## cell barcode & umi extraction-----------------------------------------------
mkdir umi_extract

# extract reads amplified from TT strategy
cat cell_list.txt | awk '$1~"TT"{print $1}' | while read sample; do

cutadapt -j 24 -m 26 -M 350 -l 28 --rc --report=minimal \
-g "ACACTCTTTCCCTACACGACGCTCTTCCGATCT;e=0.2;o=20...TTTTTTTTTTTTTTTTTTTT;e=0.2;o=20" \
-o umi_extract/${sample}.umi28.fq \
--discard-untrimmed \
raw_qual_reads/${sample}.fq

done

# extract reads amplified from TN or NN strategy
cat cell_list.txt | awk '$1~/TN|NN/{print $1}' | while read sample; do

cutadapt -j 12 -m 26 -M 30 -l 28 --rc --report=minimal \
-g "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG;e=0.2;o=20...CCTTAGCCGCTAATAGGTGAGC;e=0.2;o=20" \
-o umi_extract/${sample}.umi28.fq \
--discard-untrimmed \
raw_qual_reads/${sample}.fq

done



# umi cluster ----------------------------------------------------
mkdir umi_cluster

cat cell_list.txt | cut -f1 | while read sample; do

usearch -cluster_fast umi_extract/${sample}.umi28.fq \
-id 0.88 -mincols 24 -sizeout -sort size -strand plus \
-centroids umi_cluster/${sample}.centroids.fa \
-consout umi_cluster/${sample}.cons.fa \
-uc umi_cluster/${sample}.uc

done



# generate bc-umi table for R analysis-------------------------------
rm sample_umi_cluster_info.txt
cat cell_list.txt | cut -f1 | while read sample; do

cat umi_cluster/${sample}.centroids.fa | paste - - |
sed -e 's/^>//g' -e 's/ rc//g' -e 's/size=//g' -e 's/;/\t/g' | cut -f1,2,4 > temp

cat umi_cluster/${sample}.cons.fa | paste - - |
sed -e 's/^>//g' -e 's/size=//g' -e 's/;/\t/g' | cut -f1 |
paste - temp | awk -v var=${sample} '{print var,$0}' OFS="\t" >> sample_umi_cluster_info.txt

done



# reads bin--------------------------------------------
mkdir reads_bin; cd reads_bin

cat ../cell_list.txt | cut -f1 | while read sample; do

awk -v v=${sample} '{if($1~/[H|S]/)print">"$9"\t"v"."$2".fa"}' ../umi_cluster/${sample}.uc > ${sample}.hit
cat ../clean_qual_reads/${sample}.fa | paste - - |
awk 'NR==FNR{a[$1]=$2}NR>FNR{if($1 in a)print $0"\t"a[$1]}' ${sample}.hit - |
awk '{print $1"\n"$2 > $3}'

done


# reads-bin cluster -------------------------------------------
mkdir ../reads_cluster/

ls *fa | xargs -I {} -P 32 sh -c '
usearch -cluster_fast {} \
-id 0.99 -maxdiffs 5 -strand plus -sizeout -sort size \
-centroids ../reads_cluster/{} > /dev/null 2>&1
'


# stat of reads bin cluster ------------------------------
cd ../reads_cluster; rm ../reads_cluster.stat

ls *fa | xargs -I {} -P 100 sh -c '
echo {} $(grep -c "=" {}) $(grep "=" {} | cut -d "=" -f2 |tr "; " " ") \
>> ../reads_cluster.stat
'

# singlet reads-bin 
awk -F " " '{if($2==1 && $3==1)print$0}' reads_cluster.stat | wc -l

# potential chimera reads-bin 
awk -F " " '{if($2>1 && $3/$4<3)print$0}' reads_cluster.stat | wc -l

mkdir chimera_export
awk -F " " '{if($3/$4<3)print$0}' reads_cluster.stat | cut -d " " -f1 |
xargs -I {} cp reads_cluster/{} chimera_export/{}



# export dominate cluster ------------------------------

cat reads_cluster.stat |
awk -F " " '{if($2==1){print$1}else if($2>1 && $3/$4>=4){print$1}}' |
sed 's/.fa//' | xargs -I {} -P 100 sh -c '

echo {} $(head -1 reads_cluster/{}.fa) |
sed -e "s/[.]/ Cluster/3" -e "s/>//g" -e "s/;//g" -e "s/size=/ /g" |
sed -e "s/ /\t/g" >> dominant_cluster.txt
'



# Assigning cDNA barcode to single cell in R, generate pb_sc_match.txt file
# sc reads cluster bin & cluster-------------------------------------
mkdir sc_cluster

cat pb_sc_match.txt | sed 1d | xargs -I {} -P 100 sh -c '

info=({}); sample=${info[0]}; plant=${info[1]}; lib=${info[2]}; sc=${info[5]}; clu=${info[9]/Cluster/}
cat reads_cluster/${sample}.${clu}.fa | tr "\n" "\t" | sed "s/\t>/\n>/g" |
head -1 | tr "\t" "\n" >> sc_cluster/${plant}.${lib}.${sc}.fa
'

ls sc_cluster/*.fa | while read i; do

usearch -cluster_fast ${i} \
-id 0.99 -maxdiffs 10 -strand plus -sizein -sizeout -sort size \
-centroids ${i/fa/n2.centroid} \
-consout ${i/fa/n2.cons} > /dev/null 2>&1

done


# extract representative consensus 
ls sc_cluster/*n2.cons | while read i; do

cat ${i} | tr "\n" "\t" | sed 's/\t>/\n>/g' | head -1 |
sed -e 's/\t/\n/' -e 's/\t//g' -e 's/$/\n/' | sed '/^$/d' |
sed -e "s/>/>${i}_/g" -e 's/;//g' -e 's/size=/_N/' -e 's/n2.//' \
>> cell.n2.fa

done



# remapping to generate cstag
minimap2 -cx map-pb --cs=long -t 12 -A 2 -B 2 -O 2,12 --score-N 0 ref.mmi cell.n2.fa |
awk '{if($8==0 && $9==1030) print $0}' |
cut -f1,24 > cell.n2.cstag


# generate mutation table
cat cell.n2.cstag | cut -f2 |
sed -e 's/cs:Z://g' -e 's/=//g' | tr "ATCG" "0000" | sed 's/+[atcg]*//g' |
sed -e 's/*at/A/g' -e 's/*ac/B/g' -e 's/*ag/C/g'|
sed -e 's/*ta/D/g' -e 's/*tc/E/g' -e 's/*tg/F/g'|
sed -e 's/*ca/G/g' -e 's/*ct/H/g' -e 's/*cg/I/g'|
sed -e 's/*ga/J/g' -e 's/*gt/K/g' -e 's/*gc/L/g'|
sed -e 's/*[atcg]n/N/g' -e 's/[atcg]/d/g' -e 's/-//g' |
paste cell.n2.cstag - | cut -f1,3 |
awk 'BEGIN{printf("ref\t%-1030s\n",0)}{print $0}' | tr " " "0" |
awk '{gsub(/./,"& ",$2);print}' > cell.n2.digital.txt


# generate input file for tree
cat cell.n2.cstag| cut -f2 |
sed -e 's/cs:Z://g' -e 's/=//g' | tr "ATCG" "0000" | sed 's/+[atcg]*//g' |
sed -e 's/*ga/1/g' -e 's/*ct/1/g' |
sed 's/*[atcg][atcgn]/0/g' |sed -e 's/[atcg]/0/g' -e 's/-//g' |
paste cell.n2.cstag - | cut -f1,3 |
sed -e 's/^/>/g' | sed -e 's/\t/\n/g' |
awk 'BEGIN{printf(">ref\n%-1030s\n",0)}{print $0}' |
tr " " "0" > cell.n2.bin.fa


# iqtree
iqtree -s cell.n2.bin.fa -m GTR2+FO+G4 -B 1000 -nt 4

