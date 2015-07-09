import subprocess
import os
import click

from operator import itemgetter
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
            ['gzip', fasta_name]
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
    lines = [line.strip().split('\t') for line in contents
              if line.startswith(chrom + '\t') == True]
    indels_unsorted = [[i[0], int(i[1])] + i[2:] for i in lines]
    indels = sorted(indels_unsorted, key=itemgetter(1))
    map = []
    dodgy = []
    # start from end of chromosome to retain coordinates
    for indel in reversed(indels):
        # ensure cooordinate is within chromosome length
        if indel[1] <= len(fasta):
            indel_start = indel[1] -1
            indel_end = indel[1] + len(indel[3]) -1
            if fasta[indel_start:indel_end] == indel[3]:
                fasta = fasta[:indel_start] + indel[4] + fasta[indel_end:]
                map.append([indel[0],
                            str(indel[1]),
                            str(len(indel[4]) - len(indel[3]))
                            ]
                           )
            else:
                dodgy.append(
                'incorrect reference sequence: chr{} @{} VCF:{} fasta:{}'.format(
                        indel[0],
                        str(indel[1]),
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

def add_snps_indels(vcf_file, chromosome_constant, chrom_dir, genome_name):
    chromosomes = getattr(library_constants, chromosome_constant)
    chrom_files = os.listdir(chrom_dir)
    for chrom in chromosomes:
        chrom_file = [file_name for file_name in chrom_files\
                      if file_name.find('.{}.fa'.format(chrom)) > -1]
        if len(chrom_file) != 1:
            raise error('chromosome files are named incorrectly')
        else:
            fasta_seq = fasta_single_seq_parser(chrom_file[0], chrom_dir)
            indeled_fasta_seq = add_indels_with_mapping(fasta_seq,
                                                        chrom,
                                                        vcf_file)
        with open('{}_ref_with_snps-indels.Chromosome.{}.fasta'.format(
                    genome_name,
                    chrom
                    ), 'w') as outfile:
            outfile.write('>{}_chr{}\n{}\n'.format(
                            genome_name,
                            chrom,
                            indeled_fasta_seq
                            ))

def edit_and_combine_vcf(snp_vcf, indel_vcf, out_file):
    """
    Remove SNPs or INDELs that are nested within others, or overlap.
    :param snp_vcf:
    :param indel_vcf:
    :param out_file:
    :return:
    """
    with open(snp_vcf) as file:
        contents = file.readlines()
    vcf = [line.strip().split('\t') for line in contents
            if line.startswith('#') == False]
    with open(indel_vcf) as file:
        contents = file.readlines()
    vcf += [line.strip().split('\t') for line in contents
            if line.startswith('#') == False]
    colo = ['chromosome_name,start,end,strand,label']
    vcf_dict = {}
    for v in vcf:
        end = str(int(v[1])+len(v[3]))
        label = '|'.join(v[0:3] + [end] + [str(len(v[4]))])
        vcf_dict[label] = v
        colo.append(','.join(
            [v[0], # chromosome
             v[1], # start
             end, # end
             '1', # strand
             label, # header
             ]
        ))
    with open('colo_vcf.csv', 'w') as file:
        file.write('\n'.join(colo) + '\n')
    outfile = open('colo_out.csv', 'w')
    subprocess.check_call(
        ['colocator',
         '--search', 'colo_vcf.csv',
         '--target', 'colo_vcf.csv',
        ],
        stdout=outfile)
    outfile.close()
    with open('colo_out.csv') as file:
        contents = file.readlines()
    colocated = [i.strip().split(',') for i in contents]
    dup_dict = {}
    for colo in colocated:
        if colo[0] != colo[1]:
            if dup_dict.get(colo[0]):
                dup_dict[colo[0]].append(colo[1])
            else:
                dup_dict[colo[0]] = [colo[1]]
    nested = {}
    overlap = {}
    for dup in dup_dict:
        # for each entry, find upper level nest
        dup_list = sorted([dup.split('|')] + [i.split('|') for i in dup_dict[dup]],
                          key=itemgetter(1))
        widest = dup_list[0]
        for d in dup_list:
            if (int(d[1]) <= int(widest[1])) and (int(d[3]) > int(widest[3])):
                    widest = d
            # location is the same but larger INDEL is represented
            if (int(d[1]) == int(widest[1])) and (int(d[3]) == int(widest[3])):
                if int(d[4]) > int(widest[4]):
                    widest = d
        for d in dup_list:
            if (int(d[3]) > int(widest[3])) and (int(d[1]) < int(widest[3])):
                overlap['|'.join(widest)] = d
                overlap['|'.join(d)] = widest
            if (int(d[1]) > int(widest[1])) and (int(d[3]) <= int(widest[3])):
                nested['|'.join(d)] = widest
            elif (int(d[1]) >= int(widest[1])) and (int(d[3]) < int(widest[3])):
                nested['|'.join(d)] = widest
            elif (d[1] == widest[1]) and (d[3] == widest[3]):
                if d[2] != widest[2]:
                    nested['|'.join(d)] = widest
    no_collisions = []
    for vcf in vcf_dict:
        if overlap.get(vcf) == None:
            if nested.get(vcf) == None:
                no_collisions.append(vcf_dict[vcf])
    lines = ['\t'.join(i) for i in no_collisions]
    with open(out_file, 'w') as file:
        file.write('\n'.join(lines) + '\n')
    subprocess.call(['rm', 'colo_vcf.csv'])
    subprocess.call(['rm', 'colo_out.csv'])

@cli.command()
@click.argument('embl_gtf')
@click.argument('indel_map')
@click.argument('out_file')
def adjust_gtf_indel_map(embl_gtf, indel_map, out_file):
    with open(embl_gtf) as file:
        contents = file.readlines()
    anns = [i.strip().split('\t') for i in contents \
            if i.startswith('#') == False]
    with open(indel_map) as file:
        contents = file.readlines()
    map = [i.strip().split('\t') for i in contents]
    chromosomes = set([i[0] for i in map if i[0]])
    adj_gtf = []
    for chrom in chromosomes:
        chr_anns = [i for i in anns if i[0] == chrom]
        chr_map = [i for i in map if i[0] == chrom]
        chr_anns_unsorted = [i[:3] + [int(i[3]), int(i[4])] + i[5:] for i in chr_anns]
        chr_map_unsorted = [[i[0], int(i[1]), int(i[2])] for i in chr_map]
        chr_anns_sorted = sorted(chr_anns_unsorted, key=itemgetter(3))
        chr_map_sorted = sorted(chr_map_unsorted, key=itemgetter(1))
        m = 0 # starting map row
        start_adj = []
        for ann in chr_anns_sorted:
            for row in range(m,len(chr_map)-1):
                # if annotation starts before indel
                if m == 0 and chr_map_sorted[row][1] > ann[3]:
                    start_adj.append(ann[:3] + [str(ann[3])] + ann[4:])
                    break
                # if indel is less than annotation
                if chr_map_sorted[row][1] < ann[3]:
                    m += 1
                else:
                    start = str(ann[3] + chr_map_sorted[row-1][2])
                    start_adj.append(ann[:3] + [start] + ann[4:])
                    break
        chr_anns_end_sorted = sorted(start_adj, key=itemgetter(4))
        m = 0 # starting map row
        end_adj = []
        for ann in chr_anns_end_sorted:
            for row in range(m,len(chr_map)-1):
                # if annotation starts before indel
                if m == 0 and ann[4] <= chr_map_sorted[row][1]:
                    end_adj.append(ann[:4] + [str(ann[4])] + ann[5:])
                    break
                # if indel is lower than annotation, continue
                if chr_map_sorted[row][1] < ann[4]:
                    m = row
                else:
                    end = str(ann[4] + chr_map_sorted[row-1][2])
                    end_adj.append(ann[:4] + [end] + ann[5:])
                    break
        adj_gtf += end_adj
    lines = ['\t'.join(i) for i in adj_gtf]
    with open(out_file, 'w') as file:
        file.write('\n'.join(lines))

@cli.command()
@click.argument('vcf_file')
@click.argument('chromosome_constant')
@click.argument('chrom_dir')
@click.argument('genome_name')
def add_snps_indels_cli(vcf_file, chromosome_constant, chrom_dir, genome_name):
    add_snps_indels(vcf_file, chromosome_constant, chrom_dir, genome_name)
    subprocess.call(['gzip', '*.fasta'])

@cli.command()
@click.argument('snp_vcf')
@click.argument('indel_vcf')
@click.argument('out_file')
def edit_and_combine_vcf_cli(snp_vcf, indel_vcf, out_file):
    edit_and_combine_vcf(snp_vcf, indel_vcf, out_file)

if __name__ == '__main__':
    cli()