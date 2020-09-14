# Before doing any of this, have to generate .gff and sorted bams 

# convert *_18state_dense.bed 
# from this:
# chr10	0	3106200	18_Quies	0	.	0	3106200	255,255,255
# to this:
# chr1	19_DNase_chr1_237400_237600	19_DNase	237400	237600	0	.	.	19_DNase_chr1_237400_237600
export BAM=/home/FCAM/jvanoudenhove/ANALYSIS/CHIPSEQ/Final_Heart/align
export BIGBED=/home/FCAM/awilderman/ANALYSIS/Human_18_state_liftover/Heart
export OUT=/home/CAM/ewentworth/cotney/rawdata/chromatinsegmentations/superenhancers/human
export ROSE=/home/CAM/ewentworth/cotney/rawdata/chromatinsegmentations/superenhancers/human/rose
## 1 ## Select some of the states to create a rose_states.bed

# or 18-state
1_TssA
2_TssFlnk
3_TssFlnkU
4_TssFlnkD
5_Tx
6_TxWk
7_EnhG1
8_EnhG2
9_EnhA1
10_EnhA2
11_EnhWk
12_ZNF/Rpts
13_Het
14_TssBiv
15_EnhBiv
16_ReprPC
17_ReprPCWk
18_Quies

# select states

cat rose_states.txt | awk '{print "-e "$1}' > make_grep_list.sh
sed -i '1s/^/grep /' make_grep_list.sh 
sed -i -z 's/\n/ /g' make_grep_list.sh 

## 2 ## select files 
# use sorted beds! 
# sort (-k 1,1 -k 2,2n) if needed
export SEGS=/home/FCAM/jvanoudenhove/ANALYSIS/CHIPSEQ/Final_Heart/Human_embryonic_heart_pool1/ChromHMM/CHROMHMM_HUMAN_EMBRYONIC_HEART/18State_joint
cd $SEGS
ls -1 *.sorted.bed | sed 's/.sorted.bed//g' > $ROSE/file_list.txt
cd $ROSE
cat file_list.txt | awk '{print "grep -e 1_TssA -e 2_TssFlnk -e 3_TssFlnkU -e 4_TssFlnkD -e 7_EnhG1 -e 8_EnhG2 -e 9_EnhA1 -e 10_EnhA2 -e 11_EnhWk -e 14_TssBiv -e 15_EnhBiv $SEGS/"$1".sorted.bed > $ROSE/"$1"_sorted_rose_states.bed"}' > make_rose_state_beds.sh
#then run make_rose_state_beds.sh

## 3 ## Make .gff

ls -1 *_rose_states.bed > rose_bed_list.txt
for rosebed in `cat rose_bed_list.txt`
do 
	export NAME=`echo $rosebed | sed 's/.bed//g'`
	cat $rosebed| awk '{print $1"\t"$4"_"$1"_"$2"_"$3"\t"$4"\t"$2"\t"$3"\t"$5"\t"$6"\t.\t"$4"_"$1"_"$2"_"$3}' > $NAME.gff
done 
## 4 ## locate correct bams - if from our lab, make sure to take the bams that are normalized, sorted and indexed as for CHROMHMM
#see above, Jen already made these for the 18 state human heart. 
## 5 ## make rose.sh
##NEED TO RUN THIS WHOLE THING IN THE ROSE 
ls $BAM/*sorted.bam | egrep -v "half" > bam_list.txt
ls *rose_states.gff > gff_list.txt
egrep -hi "H3K27ac" bam_list.txt > H3K27ac_list.txt
egrep -hi "Input" bam_list.txt > Input_list.txt
paste gff_list.txt H3K27ac_list.txt Input_list.txt | sed 's/^/\/home\/CAM\/ewentworth\/cotney\/rawdata\/chromatinsegmentations\/superenhancers\/human\/rose\//'> ROSE_file_list.txt
echo -e "#/bin/bash\n#SBATCH --job-name=rose\n#SBATCH -N 1\n#SBATCH -n 2\n#SBATCH -c 1\n#SBATCH --partition=general\n#SBATCH --qos=general\n#SBATCH --mail-type=END\n#SBATCH --mem=48G\n#SBATCH --mail-user=wentworth@uchc.edu\n#SBATCH -o %j.out\n#SBATCH -e %j.err\nexport ROSE=/home/CAM/ewentworth/cotney/rawdata/chromatinsegmentations/superenhancers/human/rose\ncd /home/CAM/ewentworth/cotney/tools/young_computation-rose-feb35cb1d955\nmodule load python/2.7.14" | sed 's/\/bin/!\/bin/g' > hg19_rose.sh
cat ROSE_file_list.txt | awk '{print "python /home/CAM/ewentworth/cotney/tools/young_computation-rose-feb35cb1d955/ROSE_main.py -g hg19 -i "$1" -r "$2" -o $ROSE -c "$3}' >> hg19_rose.sh
sbatch hg19_rose.sh


## 6 ## When finished create bed for each tissue/timepoint
for sample in *SuperEnhancers.table.txt
do
export NAME=`echo $sample | sed 's/_18state_dense_sorted_rose_states_SuperEnhancers.table.txt/_heart_superenhancer.bed/g'`
awk 'NR>6{print $2 "\t" $3 "\t" $4 "\t" $1}' $sample > $NAME
done
#now lets convert them to bigbeds for loading in the ucsc browser
#this is stored in the /home/CAM/ewentworth/cotney/rawdata/chromatinsegmentations/superenhancers/human/rose directory
module load kent-tools
for sample in *heart_superenhancer.bed
do
export NAME=`echo $sample | sed 's/.bed/.bigBed/g'`
export sorted=`echo $sample | sed 's/.bed/-sorted.bed/g'`
sort -k1,1 -k2,2n $sample > $sorted
bedToBigBed $sorted /home/CAM/ewentworth/cotney/genome/hg19/hg19.chrom.sizes $NAME
scp -r $NAME /tgc/TGCore_User_Data/WebData/cotney/hubs/ChIP/ewentworth/human_heart_18state_superenhancers
done

#output superenhancer files are found in /home/CAM/ewentworth/cotney/rawdata/chromatinsegmentations/superenhancers/human
##removed the gff/ and mappedGFF/ directories and their files to save space. starting and ending files are kept in the directory, in case they are needed or if the superenhancers need to be remade. 

