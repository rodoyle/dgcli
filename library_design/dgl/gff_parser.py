import subprocess
import os
import click
import yaml

from dgl import library_constants

@click.group()
def cli():
    pass

def parse_gff(embl_gff):
    """
    parse embl gff into yaml format
    """
    with open(embl_gff) as file:
        contents = file.readlines()
    anns = [line.strip().split(
            '\t') for line in contents if line.startswith('#') == False]
    gene_dict = {}
    for ann in anns:
        info_parsed = ann[8].split(';')
        info_dict = {}
        if ann[0].find('chr') > -1:
            chr = ann[0].split('hr')[1]
        else:
            chr = ann[0]
        for info in info_parsed:
            info_dict[info.split('=')[0]] = info.split('=')[1]
        gene_dict[ann[1]] = [ann[2], # annotation type
                             chr, # chr
                             int(ann[3]), # start
                             int(ann[4]), # end
                             info_dict] # details
    return gene_dict

@cli.command()
@click.argument('embl_gff')
@click.argument('gene_list')
def parse_gff_to_colo(embl_gff, gene_list):
    """
    parse embl gff into colocator csv format
    """
    with open(embl_gff) as file:
        contents = file.readlines()
    anns = [line.strip().split(
            '\t') for line in contents if line.startswith('#') == False]
    colo_list = {}
    for ann in anns:
        if ann[0].find('chr') > -1:
            chr = ann[0].split('hr')[1]
        else:
            chr = ann[0]
        ann_info = ','.join([chr,
                             ann[3],
                             ann[4],
                             '1',
                             '_'.join([ann[1],
                                      ann[2],
                                      chr,
                                      ann[3],
                                      ann[4]])])
        if colo_list.get(chr):
            colo_list[chr].append(ann_info)
        else:
            colo_list[chr] = [ann_info]
    colo_all_list = []
    with open(gene_list) as file:
        contents = file.readlines()
    gene_info = [line.strip().split(',') for line in contents]
    chromosomes = set([gene[0] for gene in gene_info])
    for chr in chromosomes:
        with open('chr_anns.csv', 'w') as file:
            file.write('\n'.join(colo_list[chr]) + '\n')
        subprocess.call(['colocator',
                     '--search', gene_list,
                     '--target', 'chr_anns.csv',
                     '--output', 'colo_temp.csv',
                     '--ignore-strand' # only look at coordinates to judge overlap
                     ])
        with open('colo_temp.csv') as file:
            contents = file.readlines()
        colos = [line.strip().split(',') for line in contents]
        colo_all_list = colo_all_list + colos


def parse_clean_to_colo(gene_list):
    with open(gene_list) as file:
        contents = file.readlines()
    genes = [line.strip().split('\t') for line in contents]
    colo_list = []
    for gene in genes:
        colo_list.append(','.join([gene[2],
                                   gene[3],
                                   gene[4],
                                   gene[5],
                                   gene[0] + '_' + gene[1]
                                   ])
                         )
    outfile = gene_list.split('.')[0] + '_colo.csv'
    with open(outfile, 'w') as file:
        file.write('\n'.join(colo_list) + '\n')

@cli.command()
@click.argument('gene_list')
@click.argument('embl_gff')
@click.argument('outfile')
def is_ann_in_gene(gene_list, embl_gff, outfile):
    with open(gene_list) as file:
        contents = file.readlines()
    genes = [line.strip().split('\t') for line in contents]
    with open(embl_gff) as file:
        contents = file.readlines()
    anns = [line.strip().split(
            '\t') for line in contents if line.startswith('#') == False]
    chromosomes = set([gene[2] for gene in genes])
    ann_in_gene = []
    for chr in chromosomes:
        chr_anns = [i for i in anns if i[0] == 'chr' + chr]
        chr_genes = [i for i in genes if i[2] == chr]
        for gene in chr_genes:
            for ann in chr_anns:
                if gene[3] <= ann[3]:
                    if gene[4] >= ann[4]:
                        ann_in_gene.append('\t'.join(ann + gene[0:2]))
    with open(outfile, 'w') as file:
        file.write('\n'.join(ann_in_gene) + '\n')




if __name__ == '__main__':
    cli()