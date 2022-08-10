import pysam,sys,csv
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

#load reference sequence
ref = str(snakemake.input['reference'])
print(ref)

for seq_record in SeqIO.parse(ref, "fasta"):
    print(seq_record.id)
    ref_sequence = seq_record.seq
ref_aa = ref_sequence.translate()

def vcf_to_aa(ref_sequence,pos,alt):
    """
    Get amino acid information for every SNV
    """

    aa_pos = pos // 3
    codon_pos = pos % 3

    ref_aa_seq = ref_sequence.translate()
    ref_aa = ref_aa_seq[aa_pos]

    sub_ref_seq = ref_sequence[pos-codon_pos:pos-codon_pos+3]

    x = sub_ref_seq.tomutable()
    x[codon_pos] = str(alt)
    y = x.toseq()
    alt_aa = y.translate()[0]

    return(ref_aa,alt_aa)


def vcf_to_ext_df(data_path, out_path):

    myvcf = pysam.VariantFile(data_path,'r')

    pos = []
    ref = []
    alt = []
    pos_alt = []
    freq = []
    aa_pos = []
    aa_ref = []
    aa_alt = []

    for var in myvcf:
        for i in range(0,len(var.alts)):
            pos.append(var.pos - 42)
            ref.append(var.ref)
            alt.append(list(var.alts)[0])
            freq.append([var.info['AF']][i])
            pos_alt.append(str(var.pos-42) + "_" + list(var.alts)[0])
            aa = vcf_to_aa(ref_sequence=ref_sequence,pos=var.pos-1,alt=list(var.alts)[0])
            aa_pos.append((var.pos - 42 + 2) // 3)
            aa_ref.append(aa[0])
            aa_alt.append(aa[1])



    df = pd.DataFrame(
    {
        "pos":pos,
        "ref":ref,
        "alt":alt,
        "pos_alt":pos_alt,
        "freq":freq,
        "aa_pos":aa_pos,
        "aa_ref":aa_ref,
        "aa_alt":aa_alt
    })

    #classify mutations

    df["SynNonSyn"] = np.where(df["aa_ref"] == df["aa_alt"],"Syn","NonSyn")
    df["MutType"] = df[['ref','alt']].apply(lambda x: ''.join(x),axis=1)

    df.to_csv(out_path, sep=',')

vcf_to_ext_df(snakemake.input[0], snakemake.output[0])
