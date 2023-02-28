import pysam
import os
import argparse
import pandas as pd
parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bam', type=str,
                    help='bam input')
parser.add_argument('-o', '--output', type=str,
                    help='output file directory')
parser.add_argument('-t', '--tag', type=str,
                    help='BAM tag')
parser.add_argument('-r', '--reference', type=str,
                    help='Input barcode list')



if __name__ == "__main__":
    # parsing arguments
    args = parser.parse_args()
    input_file = pysam.AlignmentFile(args.bam)
    outdir = args.output;
    prefix = args.bam.split('/')[-1].split(".bam")[0]
    barcode_list = args.reference
    barcode_list = set(list(pd.read_csv(barcode_list, header=None).iloc[:,0]))
    tag = args.tag

    created_sam = set()
    opened_sam = {}
    opened_sam_queue = []
    max_queue_length = 100


    for read in input_file.fetch():
        name = int(read.qname.split(".")[-1])
        bam_fp = None
        try:
            barcode = read.get_tag(tag)
        except KeyError:
            continue
        if barcode not in barcode_list:
            continue
        if barcode not in opened_sam.keys():
            if len(opened_sam_queue) == max_queue_length:
                poped_name = opened_sam_queue.pop()
                fp = opened_sam.pop(poped_name)
                fp.close()
            if barcode not in created_sam:
                opened_sam[barcode] = pysam.AlignmentFile(os.path.join(outdir, "{}.{}.demultiplexed.sam".format(prefix, barcode)), 'w', template = input_file)
                opened_sam[barcode].close()
                opened_sam[barcode] = open(os.path.join(outdir, "{}.{}.demultiplexed.sam".format(prefix, barcode)), "a")
                created_sam.add(barcode)
            else:
                opened_sam[barcode] = open(os.path.join(outdir, "{}.{}.demultiplexed.sam".format(prefix, barcode)), "a")
            opened_sam_queue.append(barcode)

        bam_fp = opened_sam[barcode]
        bam_fp.write(read.to_string()+'\n')

    for fp in opened_sam.values():
        fp.close()
    for out_sam in os.listdir(outdir):
        pysam.sort("-o", os.path.join(outdir,out_sam.replace("sam","bam")), os.path.join(outdir,out_sam))
        os.system("rm " + os.path.join(outdir,out_sam))