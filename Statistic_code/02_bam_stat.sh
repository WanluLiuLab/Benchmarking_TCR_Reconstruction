##candidat reads statistics
filelist=$(cat file_head.txt)
for filehead in $filelist
do
	###01 candidate reads
	rm ./02_bam_statistic/${filehead}/01_totalreads_count.txt
	mkdir ./02_bam_statistic/${filehead}
	wc -l 01_bam_output/${filehead}/*toassemble.fq | awk -F " " '{print $1}' > ./02_bam_statistic/${filehead}/01_candidatereads_count.txt
	sed -i '$d' ./02_bam_statistic/${filehead}/01_candidatereads_count.txt
	wc -l 01_bam_output/${filehead}/*toassemble.fq | awk -F "." '{print $2}' > ./02_bam_statistic/${filehead}/01_cellbarcode.txt
	sed -i '$d' ./02_bam_statistic/${filehead}/01_cellbarcode.txt
	filename=$(ls ../00_data/02_bam_split/${filehead}/*.bam)
	for f in $filename
	do
		#echo $f
		count=$(samtools view -c $f)
		echo $count >> ./02_bam_statistic/${filehead}/01_totalreads_count.txt
	done
	paste ./02_bam_statistic/${filehead}/01_cellbarcode.txt ./02_bam_statistic/${filehead}/01_candidatereads_count.txt ./02_bam_statistic/${filehead}/01_totalreads_count.txt > ./02_bam_statistic/${filehead}/01_candidatereads_stat_merge.txt
	###02 assembled reads
	wc -l 01_bam_output/${filehead}/*assembled_reads.fa | awk -F " " '{print $1}' > ./02_bam_statistic/${filehead}/02_assembledreads_count.txt
	sed -i '$d' ./02_bam_statistic/${filehead}/02_assembledreads_count.txt
	wc -l 01_bam_output/${filehead}/*assembled_reads.fa | awk -F "." '{print $2}' > ./02_bam_statistic/${filehead}/02_assembledcellbarcode.txt
	sed -i '$d' ./02_bam_statistic/${filehead}/02_assembledcellbarcode.txt
	paste ./02_bam_statistic/${filehead}/02_assembledcellbarcode.txt ./02_bam_statistic/${filehead}/02_assembledreads_count.txt > ./02_bam_statistic/${filehead}/02_assembledreads_stat_merge.txt
	###03 contigreads
	rm ./02_bam_statistic/${filehead}/03_contigreads_merge.txt
	rm ./02_bam_statistic/${filehead}/04_contigreads_annot.txt
	filename=$(ls 01_bam_output/${filehead}/*_raw.out)
	for f in $filename
	do
		#echo $f
		f1=${f##*/}
		filehead2=${f1%%_raw.out}
		barcode=${filehead2%%.demultiplexed}
		barcode2=${barcode##*.}
		echo $barcode2
		cat $f | grep ">" -A 1 --no-group-separator | grep ">" -v > ./tmp/${filehead2}_contigreads.txt
		cat $f | grep ">" | awk -F " " '{print " '$barcode2' ",$1}' > ./tmp/${filehead2}_contigname.txt
		paste ./tmp/${filehead2}_contigname.txt ./tmp/${filehead2}_contigreads.txt > ./tmp/${filehead2}_contigreads_merge.txt
		cat ./tmp/${filehead2}_contigreads_merge.txt >> ./02_bam_statistic/${filehead}/03_contigreads_merge.txt
		grep ">" 01_bam_output/${filehead}/${filehead2}_annot.fa | awk -F " " '{print " '$barcode2' ",$0}' > ./tmp/${filehead2}_contigreads_annot.txt
		cat ./tmp/${filehead2}_contigreads_annot.txt >> ./02_bam_statistic/${filehead}/04_contigreads_annot.txt
	done
done
