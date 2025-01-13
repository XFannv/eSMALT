# generate bystrand ccs reads--------------------------------------------------------------

ccs -j 24 --by-strand --min-passes 3 --min-rq 0.95 --min-length 800 --max-length 1800 \
subreads.bam bystrand_ccs.bam 

samtools fastq bystrand_ccs.bam  > ccs.fq



## mapping--------------------------------------------------------------

minimap2 -x map-pb -d ref.mmi ref.fa

mkdir mapping

cat sample_list.txt | while read sample; do

minimap2 -cx map-pb --cs=long -t 12 -A 2 -B 2 -O 2,12 ref.mmi ${sample}.fq |
awk '{if($12==60 && $10>=900 && $8==0 && $9==1030) print $0}' |
cut -f1,24 > mapping/${sample}.cstag

done



#### UMI consensus~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# extract quality reads based on mapping results-------------------------------
mkdir raw_qual_reads

cat sample_list.txt | while read sample; do

cat mapping/${sample}.cstag | awk '{print"@"$1}' > raw_qual_reads/${sample}.id

cat ${sample}.fq | paste - - - - | sed 's/ rc//g' |
awk 'NR==FNR{a[$1]=$1}NR>FNR{if($1 in a)print$0}' raw_qual_reads/${sample}.id - |
tr "\t" "\n" > raw_qual_reads/${sample}.fq

done

rm raw_qual_reads/*id


# generate quality barcode sequence based on cstag in mapping results------------------------------
mkdir clean_qual_reads

cat sample_list.txt | while read sample; do

cat mapping/${sample}.cstag | cut -f2 |
sed -e 's/cs:Z://g' -e 's/=//g' |
sed -e 's/+[atcg]*//g' -e 's/-//g' |
sed -e 's/*ga/a/g' -e 's/*ct/t/g' |
sed -e 's/*a[tcg]/a/g' -e 's/*t[acg]/t/g' -e 's/*c[ag]/c/g' -e 's/*g[tc]/g/g' |
tr atcg ATCG | paste mapping/${sample}.cstag - | cut -f1,3 |
sed -e 's/^m/>m/g' | sed -e 's/\t/\n/g' > clean_qual_reads/${sample}.fa

done


## umi extraction-------------------------------------------------
mkdir umi_extract

cat sample_list.txt | while read sample; do

cutadapt -j 24 -m 16 -M 20 --report=minimal \
-g "XCATACGAGAT;e=2;o=8...TCCCACGGCGACCACTTCG;e=4;o=18" \
-o umi_extract/${sample}.UMI.LB.fq \
--discard-untrimmed \
raw_qual_reads/${sample}.fq

done

cat sample_list.txt | while read sample; do

cutadapt -j 24 -m 16 -M 20 --report=minimal \
-g "GCTCACCTATTAGCGGCTAAGG;e=4;o=18...GATCTCGGTGX;e=2;o=8" \
-o umi_extract/${sample}.UMI.RB.fq \
--discard-untrimmed \
raw_qual_reads/${sample}.fq

done

cat sample_list.txt | while read sample; do

cat umi_extract/${sample}.UMI.LB.fq | paste - - - - > LB.temp
cat umi_extract/${sample}.UMI.RB.fq | paste - - - - > RB.temp

awk 'NR==FNR{A[$1]=$0}NR>FNR{if($1 in A){print A[$1]"\t"$0}}' LB.temp RB.temp |
awk '{print $1"\n"$2$6"\n"$3"\n"$4$8}' > umi_extract/${sample}.UMI.fq

done


# UMI cluster------------------------------------------
mkdir umi_cluster

cat sample_list.txt | while read sample; do

usearch -fastx_uniques umi_extract/${sample}.UMI.fq \
-sizeout -strand plus \
-fastaout umi_cluster/${sample}.uniq.fa

usearch -cluster_fast umi_cluster/${sample}.uniq.fa \
-id 0.8 -maxdiffs 6 -sizein -sizeout -sort size -strand plus \
-centroids umi_cluster/${sample}.centroids.fa \
-consout umi_cluster/${sample}.cons.fa \
-uc umi_cluster/${sample}.uc

done



# reads bin --------------------------------------------
mkdir reads_bin; cd reads_bin

cat ../sample_list.txt | while read sample; do

awk -v v=${sample} '{if($1~/[H|S]/)print">"$9v"."$2".fa"}' ../umi_cluster/${sample}.uc |
cut -d ";" -f1,3 | sed 's/;/\t/' > ${sample}.hit
cat ../clean_qual_reads/${sample}.fa | paste - - |
awk 'NR==FNR{a[$1]=$2}NR>FNR{if($1 in a)print $0"\t"a[$1]}' ${sample}.hit - |
awk '{print $1"\n"$2 > $3}'

done


# reads-bin cluster -------------------------------------------
mkdir ../reads_cluster/

cat ../sample_list.txt | while read sample; do

ls ${sample}*fa | xargs -I {} -P 32 sh -c '
usearch -cluster_fast {} \
-id 0.9 -maxdiffs 10 -strand plus -sizeout -sort size \
-centroids ../reads_cluster/{} > /dev/null 2>&1
'
done

# stat of reads bin cluster
cd ../reads_cluster
rm ../reads_cluster.stat

cat ../sample_list.txt | while read sample; do

ls ${sample}*fa | xargs -I {} -P 100 sh -c '
echo {} $(grep -c "=" {}) $(grep "=" {} | cut -d "=" -f2 |tr "; " " ") \
>> ../reads_cluster.stat
'
done



# export representative read of quality umi-cluster------------------------------------

cat sample_list.txt | while read sample; do

cat reads_cluster.stat | grep ${sample} | 
awk -F " " '{if($2==1 && $3>1){print$1}else if($2>1 && $3/$4>=2){print$1}}' | xargs -I {} -P 100 sh -c '
cat reads_cluster/{} | head -1 | cut -d ";" -f1 | sed "s/>//" >> temp
'
awk 'NR==FNR{a[$1]=$1}NR>FNR{if($1 in a)print$0}' temp mapping/${sample}.cstag \
> mapping/${sample}.cons.cstag

done


## tree reconstruction based on barcoding mutations----------------------------

# sampling 200 reads from each shoot branch

cat sample_list.txt | while read sample; do

cat mapping/${sample}.cons.cstag | shuf -n 200 |
sed "s/^m/${sample}_m/" >> plant.cstag

done

# generate mutation table

cat plant.cstag | cut -f2 |
sed -e 's/cs:Z://g' -e 's/=//g' | tr "ATCG" "0000" |
sed 's/+[atcg]*//g' |
sed -e 's/*at/A/g' -e 's/*ac/B/g' -e 's/*ag/C/g'|
sed -e 's/*ta/D/g' -e 's/*tc/E/g' -e 's/*tg/F/g'|
sed -e 's/*ca/G/g' -e 's/*ct/H/g' -e 's/*cg/I/g'|
sed -e 's/*ga/J/g' -e 's/*gt/K/g' -e 's/*gc/L/g'|
sed -e 's/[atcg]/d/g' -e 's/-//g' |
paste plant.cstag - | cut -f1,3 |
awk 'BEGIN{printf("ref\t%-1030s\n",0)}{print $0}' | tr " " "0" |
awk '{gsub(/./,"& ",$2);print}' > plant.digital.txt


# generate input file for iqtree

cat plant.cstag | cut -f2 |
sed -e 's/cs:Z://g' -e 's/=//g' | tr "ATCG" "0000" |
sed 's/+[atcg]*//g' |
sed -e 's/*ga/1/g' -e 's/*ct/1/g' |
sed 's/*[atcg][atcg]/0/g' |
sed -e 's/[atcg]/0/g' -e 's/-//g' |
paste plant.cstag - | cut -f1,3 |
sed -e 's/^/>/g' | sed -e 's/\t/\n/g' |
awk 'BEGIN{printf(">ref\n%-1030s\n",0)}{print $0}' |
tr " " "0" > plant.bin.fa


# iqtree
iqtree -s plant.bin.fa -m GTR2+FO+I+G4 -B 1000 -nt 4

