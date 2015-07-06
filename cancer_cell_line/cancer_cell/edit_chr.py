import subprocess
import os
import click

from dgl import library_constants

@click.group()
def cli():
    pass

def fasta_single_seq_parser(fasta_file, fasta_dir):
    """
    return a single fasta sequence as a string for a single-sequence fasta file
    e.g. a chromosome sequence
    :param fasta_file: fasta file name
    :return:
    """
    if fasta_file.endswith('.gz'):
        zipped = True
        subprocess.call(
            ['gzip', '-d', fasta_dir + fasta_file]
        )
        fasta_name = fasta_dir + fasta_file.split('.gz')[0]
    else:
        zipped = False
        fasta_name = fasta_dir + fasta_file
    with open(fasta_name) as file:
        fasta_lines = file.readlines()
    if zipped:
        subprocess.call(
            ['gzip', fasta_dir + fasta_name]
        )
    fasta_seq = ''
    for line in fasta_lines[1:]:
        fasta_seq += line.strip()
    return fasta_seq

def add_snps(fasta, chrom, snp_vcf):
    """
    fasta is string input of chromsome seqence
    chrom is string of chromosome name
    snp_vcf is vcf file name for snps
    use caveman.vcf-type file to apply SNPs to
    ref chromosome sequence
    """
    dodgy = []
    with open(snp_vcf) as file:
        contents = file.readlines()
    snps = [line.strip().split('\t') for line in contents
            if line.startswith(chrom + '\t') == True]
    for snp in snps:
        if int(snp[1]) <= len(fasta):
            if fasta[int(snp[1])-1] == snp[3]:
                fasta = fasta[0:int(snp[1])-1] + snp[4] + fasta[int(snp[1]):]
            else:
                dodgy.append(
                    'incorrect reference sequence: chr{} @{} VCF:{} fasta:{}'.format(
                        snp[0],
                        snp[1],
                        snp[2],
                        fasta[int(snp[1])-1]))
    if len(dodgy) > 0:
        with open('incorrect_snps.txt', 'a') as error_file:
            error_file.write('\n'.join(dodgy) + '\n')
    return fasta

def add_indels_with_mapping(fasta, chrom, indel_vcf):
    """
    fasta is string input of chromsome seqence
    chrom is string of chromosome name
    indel_vcf is vcf file name for indels
    use pindel.vcf-type file to integrate INDELs into
    chromosome sequence and make a  cumulative map file
    to adjust feature coordinates.
    """
    with open(indel_vcf) as file:
        contents = file.readlines()
    indels = [line.strip().split('\t') for line in contents
              if line.startswith(chrom + '\t') == True]
    map = []
    dodgy = []
    # start from end of chromosome to retain coordinates
    for indel in reversed(indels):
        # ensure cooordinate is within chromosome length
        if int(indel[1]) <= len(fasta):
            indel_start = int(indel[1]) -1
            indel_end = int(indel[1]) + len(indel[3]) -1
            if fasta[indel_start:indel_end] == indel[3]:
                fasta = fasta[:indel_start] + indel[4] + fasta[indel_end:]
                map.append([indel[0],
                            indel[1],
                            str(len(indel[4]) - len(indel[3]))
                            ]
                           )
            else:
                dodgy.append(
                'incorrect reference sequence: chr{} @{} VCF:{} fasta:{}'.format(
                        indel[0],
                        indel[1],
                        indel[3],
                        fasta[indel_start:indel_end]))
    if len(dodgy) > 0:
        with open('incorrect_indels.txt', 'a') as error_file:
            error_file.write('\n'.join(dodgy) + '\n')
    map_cum = []
    for change in reversed(map):
        if len(map_cum) == 0:
            map_cum.append(change)
        else:
            map_cum.append([change[0],
                            change[1],
                            str(int(map_cum[len(map_cum)-1][2]) + int(change[2]))
                            ]
                           )
    map_line = ['\t'.join(entry) for entry in map_cum]
    with open('indel_map.txt', 'a') as map_file:
        map_file.write('\n'.join(map_line) + '\n')
    return fasta

def add_snps_indels(snp_vcf, indel_vcf, chromosome_constant, chrom_dir, genome_name):
    chromosomes = getattr(library_constants, chromosome_constant)
    chrom_files = os.listdir(chrom_dir)
    for chrom in chromosomes:
        chrom_file = [file_name for file_name in chrom_files\
                      if file_name.find('.{}.fa'.format(chrom)) > -1]
        if len(chrom_file) != 1:
            raise error('chromosome files are named incorrectly')
        else:
            fasta_seq = fasta_single_seq_parser(chrom_file[0],
                                                chrom_dir)
            snped_fasta_seq = add_snps(fasta_seq,
                                       chrom,
                                       snp_vcf)
            indeled_fasta_seq = add_indels_with_mapping(snped_fasta_seq,
                                                        chrom,
                                                        indel_vcf)
        with open('{}_ref_with_snps-indels.Chromosome.{}.fasta'.format(
                    genome_name,
                    chrom
                    ), 'w') as outfile:
            outfile.write('>{}_chr{}\n{}\n'.format(
                            genome_name,
                            chrom,
                            indeled_fasta_seq
                            ))

@cli.command()
@click.argument('snp_vcf')
@click.argument('indel_vcf')
@click.argument('chromosome_constant')
@click.argument('chrom_dir')
@click.argument('genome_name')
def add_snps_indels_cli(snp_vcf, indel_vcf, chromosome_constant, chrom_dir, genome_name):
    add_snps_indels(snp_vcf, indel_vcf, chromosome_constant, chrom_dir, genome_name)
    subprocess.call('gzip', '*.fasta')

if __name__ == '__main__':
    cli()