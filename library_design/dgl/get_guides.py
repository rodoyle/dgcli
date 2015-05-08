import os
import subprocess
import re
import math
import click

@click.group()
def cli():
    pass

@cli.command()
@click.argument('gene_details')
@click.argument('biodata_prefix')
@click.argument('outfile')
@click.argument('chromosome-list', 'chromosome_constant')
def get_chr_guides(gene_details, biodata_prefix, outfile, chromosome_constant):
    """
    find all guides in gene list from
    chromosome fasta sequences
    """
    chromosomes = getattr(library_constants, chromosome_constant)
    guide_pattern = getattr(library_constants, GUIDE_PATTERN)
    reverse_pattern = getattr(library_constants, REVERSE_GUIDE_PATTERN)
    with open(gene_details) as file:
        contents = file.readlines()
    all_genes = [line.strip().split('\t') for line in contents]
    guide_list = []
    for chr in chromosomes:
        chr_genes = [gene for gene in all_genes if gene[1] == chr]
        if chr_genes:
            with open(biodata_prefix + chr + '.fasta') as file:
                sequence = file.read()
            chr_fasta = sequence.strip().split('\n')[1]
            del sequence
            for gene in chr_genes:
                gene_seq = chr_fasta[int(gene[2])-1:int(gene[3])]
                fwd_guide_seqs = re.findall(guide_pattern, gene_seq)
                rev_guide_seqs = re.findall(reverse_pattern, gene_seq)
                fwd_guide_loci = re.finditer(guide_pattern, gene_seq)
                rev_guide_loci = re.finditer(reverse_pattern, gene_seq)
                for n,position in enumerate(fwd_guide_loci):
                    guide_list.append([gene[0], gene[1],
                                       str(int(gene[2])+position.start(1)+1),
                                       '1',
                                       fwd_guide_seqs[n]])
                for n,position in enumerate(rev_guide_loci):
                    guide_list.append([gene[0], gene[1],
                                       str(int(gene[2])+position.start(1)+1),
                                       '-1',
                                       rev_guide_seqs[n]])
    lines = []
    for i in guide_list:
        lines.append('\t'.join(i))
    with open(outfile, 'w') as file:
        file.write('\n'.join(lines))

@cli.command()
@click.argument('cds_list')
@click.argument('guide_list')
@click.argument('gene_list')
@click.argument('outfile')
def is_coding(cds_list, guide_list, gene_list, outfile):
    """
    determine if guides cut in CDS regions and outputs
    guides in coding region only.
    also adds unique guide ID
    """
    with open(cds_list) as file:
        contents = file.readlines()
    all_cdss = [line.strip().split('\t') for line in contents]
    with open(guide_list) as file:
        contents = file.readlines()
    all_guides = [line.strip().split('\t') for line in contents]
    with open(gene_list) as file:
        contents = file.readlines()
    genes = [line.strip().split('\t') for line in contents]
    coding_guides = []
    columns = len(all_guides[0])
    for gene in genes:
        gene_guides = [guide for guide in all_guides if guide[0] == gene[0]]
        gene_cdss = [cds for cds in all_cdss if cds[0] == gene[0]]
        for guide in gene_guides:
            if guide[3] == '1':
                cut_site = int(guide[2]) + 21
            elif guide[3] == '-1':
                cut_site = int(guide[2]) + 9
            for cds in gene_cdss:
                if len(guide) == columns:
                    if cut_site > int(cds[4]) and cut_site <= int(cds[5]):
                        guide.append('coding')
                        coding_guides.append(guide[0:columns])
    lines = []
    for i in coding_guides:
        lines.append('_'.join(i[0:columns-1]) + '\t' + '\t'.join(i))
    with open(outfile, 'w') as file:
        file.write('\n'.join(lines))

if __name__ == '__main__':
    cli()