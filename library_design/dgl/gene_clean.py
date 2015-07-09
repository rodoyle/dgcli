import click

STRAND = {'+':'1',
          '-':'-1'}

@click.group()
def cli():
    pass

def parse_gtf(embl_gtf):
    """
    parse embl gtf into yaml format
    keys:
    gene = gene_id
    transcript = transcript_id
    exon = transcript_id|exon_id
    CDS = transcript_id|protein_id|exon_number
    """
    with open(embl_gtf) as file:
        contents = file.readlines()
    anns = [line.strip().split(
            '\t') for line in contents if line.startswith('#') == False]
    strand = {'+':'1',
              '-':'-1'}
    gene_dict = {}
    for ann in anns:
        info_parsed = ann[8].split('; ')
        info_dict = {}
        for info in info_parsed:
            info_dict[info.split(' "')[0]] = info.split('"')[1]
        ann_info = [ann[2], # annotation type
                    ann[0], # chr
                    int(ann[3]), # start
                    int(ann[4]), # end
                    strand[ann[6]], # strand
                    info_dict] # details
        if ann[2] == 'gene':
            gene_dict[info_dict['gene_id']] = ann_info
        elif ann[2] == 'transcript':
            gene_dict[
                '{}|{}'.format(info_dict['gene_id'],
                               info_dict['transcript_id'],
                                )
                     ] = ann_info
        elif ann[2] == 'exon':
            gene_dict['{}|{}'.format(info_dict['transcript_id'],
                                     info_dict['exon_id'],
                                     )
                     ] = ann_info
        elif ann[2] == 'CDS':
            gene_dict['{}|{}|{}'.format(info_dict['transcript_id'],
                                        info_dict['protein_id'],
                                        info_dict['exon_number'],
                                        )
                     ] = ann_info
    return gene_dict

def get_annotations_gtf(embl_gtf, gene_list, embl):
    """
    Extract gene and CDS annotations for all genes
    in the gene_list from embl .gtf file.
    Gene list is a tab-delimited file from Excel with one column
    embl = are gene names embl ids? - if False, uses HUGO names as input
    """
    ann_dict = parse_gtf(embl_gtf)
    with open(gene_list) as file:
        contents = file.read()
    genes = set(contents.strip().split('\n'))
    if len(genes) == 1:
        genes = set(contents.strip().split('\r'))
    cds_list = []
    anno_list = []
    not_found = []
    for gene in genes:
        gene_details = []
        for anno in ann_dict:
            if len(ann_dict[anno][1]) < 3: # on chr not non-chr contigs
                if embl == False:
                    if ann_dict[anno][0] == 'gene':
                        if ann_dict[anno][5].get('gene_name'):
                            if ann_dict[anno][5]['gene_name'] == gene:
                                gene_details.append([gene,
                                                     ann_dict[anno][5]['gene_id'],
                                                     ann_dict[anno][1],
                                                     str(ann_dict[anno][2]),
                                                     str(ann_dict[anno][3]),
                                                     ann_dict[anno][4]])
                    elif ann_dict[anno][0] == 'CDS':
                        if ann_dict[anno][5].get('gene_name'):
                            if ann_dict[anno][5]['gene_name'] == gene:
                                cds_list.append([gene,
                                                 anno,
                                                 ann_dict[anno][1],
                                                 str(ann_dict[anno][2]),
                                                 str(ann_dict[anno][3]),
                                                 ann_dict[anno][4]])
                else:
                    if ann_dict[anno][0] == 'gene':
                        if ann_dict[anno][5]['gene_id'] == gene:
                            gene_details.append([ann_dict[anno][5]['gene_name'],
                                                 gene,
                                                 ann_dict[anno][1],
                                                 str(ann_dict[anno][2]),
                                                 str(ann_dict[anno][3]),
                                                 ann_dict[anno][4]])
                    elif anno[2] == 'CDS':
                        if ann_dict[anno][5]['gene_id'] == gene:
                            cds_list.append([gene,
                                             anno,
                                             ann_dict[anno][1],
                                             str(ann_dict[anno][2]),
                                             str(ann_dict[anno][3]),
                                             ann_dict[anno][4]])
        if len(gene_details) == 0:
            not_found.append(gene)
        else:
            anno_list += gene_details
    lines = []
    for i in anno_list:
        lines.append('\t'.join(i))
    file_name = gene_list.split('.')[0] + '_clean.txt'
    with open(file_name, 'w') as file:
        file.write('\n'.join(lines))
    lines = []
    for i in cds_list:
        lines.append('\t'.join(i))
    file_name = gene_list.split('.')[0] + '_CDSs.txt'
    with open(file_name, 'w') as file:
        file.write('\n'.join(lines))
    file_name = gene_list.split('.')[0] + '_not_found.txt'
    with open(file_name, 'w') as file:
        if len(not_found) == 0:
            file.write('No Missing Genes')
        else:
            file.write('\n'.join(not_found))

@cli.command()
@click.argument('embl_gtf')
@click.argument('gene_list')
@click.argument('embl', default=False)
def get_annotations_gtf_cli(embl_gtf, gene_list, embl):
    get_annotations_gtf(embl_gtf, gene_list, embl)

def find_ensembl_from_alias_human(gene_list):
    with open('/home/neil/biodata/human_protein_coding_gene_alias.txt') as file:
        contents = file.readlines()
    alias = [i.strip().split('\t') for i in contents]
    with open(gene_list) as file:
        contents = file.read()
    genes = contents.strip().split('\n')
    ensembl = []
    multiple = []
    not_found = []
    for gene in genes:
        ensg = []
        for ale in alias[1:]:
            if gene in ale:
                if ale[19].startswith('ENSG'):
                    ensg.append([gene, ale[19]])
                else:
                    not_found.append(gene)
                if len(ensg) != 0:
                    multiple.append(gene)
        if len(ensg) == 0:
            not_found.append(gene)
        else:
            ensembl += ensg
    file_name = gene_list.split('.')[0] + '_multiple_matches.txt'
    with open(file_name, 'w') as file:
        file.write('\n'.join(multiple) + '\n')
    file_name = gene_list.split('.')[0] + '_not_found.txt'
    with open(file_name, 'w') as file:
        file.write('\n'.join(not_found) + '\n')
    lines = '\n'.join([i[1] for i in ensembl])
    file_name = gene_list.split('.')[0] + '_ensg_ids.txt'
    with open(file_name, 'w') as file:
        file.write(lines + '\n')
    embl_gtf = '/home/neil/biodata/embl_annotations/Homo_sapiens.GRCh38.79.gtf'
    get_annotations_gtf(embl_gtf, file_name, embl=True)

@cli.command()
@click.argument('gene_list')
def find_ensembl_from_alias_human_cli(gene_list):
    find_ensembl_from_alias_human(gene_list)






if __name__ == '__main__':
    cli()