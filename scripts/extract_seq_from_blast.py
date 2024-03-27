#!/usr/bin/env python3

import os
import subprocess
from io import StringIO

import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def bedtools_cmd(sn, fn, seq, start, end, frame):
    """Calls bedtools getfasta to extract the best best blast hit from the input genome assemblies"""
    strand = {True: '+', False: '-'}[int(frame) > 0]
    start, end = {'-': (end,start), '+': (start, end)}[strand]
    start = str(int(start)-1)
    bedline = [seq, start, end, ".", frame, strand]
    bedline = '\t'.join(bedline)
    with open(sn+".bed", 'w') as bedF:
        bedF.write(bedline+'\n')

    cmdline = ["bedtools", "getfasta", "-fi", fn, "-bed", sn+".bed", "-s"]

    result = subprocess.run(cmdline, stdout=subprocess.PIPE).stdout.decode('utf-8')
    os.remove(sn+".bed")

    return result

def process_blast_result(samplename, sampledict, header_dict,idx=0):
    """Parse tblastn output table for a sample and determine the "best" (longest) hit
    according to the following rule:
        - Sample is not a genome assembly: get the overall best hit
        - Sample is a genome assembly:
            - Has a clearly labeled mitochondrial contig: get the best hit
              among the hits against the mitochondrial contig
            - else: get the overall best hit
    """
    sample_df = pd.read_csv(samplename+"_tblastn.tsv", sep='\t', comment='#', header=None)

    sample = sampledict[samplename]
    mode = sample['mode']

    if mode != 'genome':
        best_hit = sample_df.sort_values(8, ascending=False).iloc[idx,:].values.tolist()
        extracted = bedtools_cmd(samplename, sample['path'], str(best_hit[0]), str(best_hit[2]), str(best_hit[3]), str(best_hit[4]))
        return (sample_df, extracted)

    hdrs = header_dict[samplename]
    mito = [x for x in hdrs if x.find('mito') != -1]

    if len(mito) < 1:
        best_hit = sample_df.sort_values(8, ascending=False).iloc[idx,:].values.tolist()
        extracted = bedtools_cmd(samplename, sample['path'], str(best_hit[0]), str(best_hit[2]), str(best_hit[3]), str(best_hit[4]))
        return (sample_df, extracted)

    mitocontig = mito[0].split()[0]
    mito_df = sample_df[sample_df[0] == mitocontig]
    best_hit = mito_df.sort_values(8, ascending=False).iloc[idx,:].values.tolist()
    extracted = bedtools_cmd(samplename, sample['path'], str(best_hit[0]), str(best_hit[2]), str(best_hit[3]), str(best_hit[4]))

    return (sample_df, extracted)

def write_fasta_from_blast(samplename, sampledict, header_dict):
    df, blast_best_hit = process_blast_result(samplename, sampledict, header_dict)
    fasta_io = StringIO(blast_best_hit)
    this_seq = SeqIO.read(fasta_io, "fasta")
    this_seq.id=samplename
    this_seq_prot = str(this_seq.seq.translate())
    this_seq_prot_mito = str(this_seq.seq.translate(table='5'))
    if this_seq_prot_mito.find("*") != -1:
        print(f"Extracted CDS has a stop codon!!! for sample {samplename}")
        return 1

    SeqIO.write(this_seq, f"{samplename}_COI_nt.fasta", "fasta")

    this_seq.seq = this_seq.seq.translate(table='5')
    SeqIO.write(this_seq, f"{samplename}_COI_aa.fasta", "fasta")

    return 0
    
def main():
    samples = {}

    p = "genomes.txt"
    if os.path.exists(p):
        with open(p) as genomeF:
            genome_samples = [x.strip() for x in genomeF.readlines()]
        for sample in genome_samples:
            k = sample.split("/")[-1]
            k = k.replace(".fa", "")
            k = k.replace(".fna", "")
            thisElement = {'path': sample, 'mode': 'genome'}
            samples[k] = thisElement

    p = "transcriptomes.txt"
    if os.path.exists(p):
        with open(p) as transF:
            transcriptome_samples = [x.strip() for x in transF.readlines()]
        for sample in transcriptome_samples:
            k = sample.split("/")[-1]
            k = k.replace(".fa", "")
            k = k.replace(".fna", "")
            thisElement = {'path': sample, 'mode': 'transcriptome'}
            samples[k] = thisElement


    print("parsing genome headers...")
    genome_headers={}
    for sample in samples:
        if samples[sample]['mode'] != 'genome':
            continue
        genomeF = samples[sample]['path']
        headers = []
        with open(genomeF, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                headers.append(record.description)
        genome_headers[sample] = headers
    print("Done.")

    for sample in samples:
        write_fasta_from_blast(sample, samples, genome_headers)

    return 0

if __name__ == "__main__":
    main()
