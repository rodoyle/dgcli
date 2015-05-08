import click

STRAND = {'+':'1',
          '-':'-1'}

@click.group()
def cli():
    pass

@cli.command()
@click.argument('embl_gtf')
@click.argument('gene_list')
@click.argument('embl', default='False')
def get_annotations_gtf(embl_gtf, gene_list, embl):
    """
    Extract gene and CDS annotations for all genes
    in the gene_list from embl .gtf file.
    Gene list is a tab-delimited file from Excel with one column
    embl = are gene names embl ids? - if False, uses HUGO names as input
    """
    with open(embl_gtf) as file:
        contents = file.readlines()
    anns = [line.strip().split('\t') for line in contents]
    with open(gene_list) as file:
        contents = file.read()
    genes = contents.strip().split('\r')
    cds_list = []
    gene_details = []
    #conditions = (
    #    anno[2] == 'gene',
    #    )
    for gene in genes:
        for anno in anns:
            if len(anno) == 9:
                if embl == False:
                    if anno[2] == 'gene':
                        g_name = anno[8].split('gene_name "')
                        if len(g_name) == 2: 
                            if g_name[1].startswith(gene + '"'):
                                gene_id = anno[8].split('"')
                                gene_details.append([gene,
                                                    gene_id[1],
                                                    anno[0],
                                                    anno[3],
                                                    anno[4],
                                                    STRAND[anno[6]]])
                    elif anno[2] == 'CDS':
                        g_name = anno[8].split('gene_name "')
                        if len(g_name) == 2: 
                            if g_name[1].startswith(gene + '"'):
                                cds_list.append(gene + '\t' + '\t'.join(anno))
                else:
                    if anno[2] == 'gene':
                        g_name = anno[8].split('ene_id "')
                        if len(g_name) == 2: 
                            if g_name[1].startswith(gene + '"'):
                                gene_details.append([gene,
                                                    anno[0],
                                                    anno[3],
                                                    anno[4],
                                                    STRAND[anno[6]]])
                    elif anno[2] == 'CDS':
                        g_name = anno[8].split('ene_id "')
                        if len(g_name) == 2: 
                            if g_name[1].startswith(gene + '"'):
                                cds_list.append(gene + '\t' + '\t'.join(anno))
    lines = []
    for i in gene_details:
        lines.append('\t'.join(i))
    file_name = gene_list.split('.')[0] + '_clean.txt'
    with open(file_name, 'w') as file:
        file.write('\n'.join(lines))
    file_name = gene_list.split('.')[0] + '_CDSs.txt'
    with open(file_name, 'w') as file:
        file.write('\n'.join(cds_list))

if __name__ == '__main__':
    cli()