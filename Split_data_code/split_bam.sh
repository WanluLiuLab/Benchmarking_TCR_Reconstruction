mkdir ./02_bam_split/CPIc_C1/
python demultiplex_bam.py -b ./01_bam_merge/CPIc_C1.bam -o ./02_bam_split/CPIc_C1/ -t CB -r ./03_barcode/CPIc_C1_barcodes.tsv
mkdir ./02_bam_split/CPIc_C2/
python demultiplex_bam.py -b ./01_bam_merge/CPIc_C2.bam -o ./02_bam_split/CPIc_C2/ -t CB -r ./03_barcode/CPIc_C2_barcodes.tsv
mkdir ./02_bam_split/CPIc_C3/
python demultiplex_bam.py -b ./01_bam_merge/CPIc_C3.bam -o ./02_bam_split/CPIc_C3/ -t CB -r ./03_barcode/CPIc_C3_barcodes.tsv
mkdir ./02_bam_split/CPIc_C4/
python demultiplex_bam.py -b ./01_bam_merge/CPIc_C4.bam -o ./02_bam_split/CPIc_C4/ -t CB -r ./03_barcode/CPIc_C4_barcodes.tsv
mkdir ./02_bam_split/CPIc_C5/
python demultiplex_bam.py -b ./01_bam_merge/CPIc_C5.bam -o ./02_bam_split/CPIc_C5/ -t CB -r ./03_barcode/CPIc_C5_barcodes.tsv
mkdir ./02_bam_split/CPIc_C6/
python demultiplex_bam.py -b ./01_bam_merge/CPIc_C6.bam -o ./02_bam_split/CPIc_C6/ -t CB -r ./03_barcode/CPIc_C6_barcodes.tsv
mkdir ./02_bam_split/CPIc_C7/
python demultiplex_bam.py -b ./01_bam_merge/CPIc_C7.bam -o ./02_bam_split/CPIc_C7/ -t CB -r ./03_barcode/CPIc_C7_barcodes.tsv

