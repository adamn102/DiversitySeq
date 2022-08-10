import pysam,sys,csv

#convert vcf data into a csv data table
#
def vcf_to_freq(data_path, out_path):

    vcf_file = data_path
    freq_out = out_path

    myvcf = pysam.VariantFile(vcf_file,'r')

    col_names = ['reference_sequence','pos','id','ref','alt','freq']

    with open(freq_out,'w') as out_file:

        writer = csv.writer(out_file)
        writer.writerow((col_names))

        for var in myvcf:
            row_data  = [var.chrom, str(var.pos), str(var.id), var.ref, str(var.alts), str(var.info['AF'])]
            writer.writerow(row_data)

vcf_to_freq(snakemake.input[0], snakemake.output[0])
