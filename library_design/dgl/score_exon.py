import click

STRAND = {'+':'1',
          '-':'-1'}

@click.group()
def cli():
    pass

@cli.command()
@click.argument('gene_list')
@click.argument('cds_list')
@click.argument('outfile')
def score_exons(gene_list, cds_list, outfile):
    """
    Output list of exons or sub-exonic regions
    within CDS with the number of CDSs they appear in.
    Input list of genes, and list of CDSs that start each line
    with the corresponding gene name.
    """
    with open(cds_list) as file:
        all_cdss = file.readlines()
    with open(gene_list) as file:
        contents = file.readlines()
    gene_list = [line.strip().split('\t') for line in contents]
    common = []
    for gene in gene_list:
        gene_cdss = [line.strip().split('\t') for line in all_cdss if line.startswith(gene[0]) == True]
        exons = []
        if gene_cdss:
            for cds in gene_cdss:
                for position in cds[4:6]:
                    exons.append(int(position))
            scores = [0] * (max(exons)-min(exons)+1)
            for exon in range(0,len(exons),2):
                for score in range(exons[exon]-min(exons), exons[exon+1]-min(exons), 1):
                    scores[score] = scores[score]+1
            max_score = max(scores)
            start = 0
            start_score = -1
            for n,score in enumerate(scores):
                if start == 0:
                    if score > 0:
                        start = n+1
                        start_score = score
                elif start > 0:
                    if score == start_score:
                        continue
                    else:
                        exon_score = float(start_score)/float(max_score)
                        common.append([gene[0], 
                                       gene[1], 
                                       str(start+min(exons)), 
                                       str(n+1+min(exons)),
                                       gene[4], 
                                       str(start_score),
                                       "%.2f" % exon_score])
                        start = 0
                        start_score = -1
    lines = []
    for i in common:
        lines.append('\t'.join(i))
    with open(outfile, 'w') as file:
        file.write('\n'.join(lines))

@cli.command()
@click.argument('exon_list')
@click.argument('guide_list')
@click.argument('outfile')
def add_exon_score(exon_list, guide_list, outfile):
    """
    add exon score column to guide list
    """
    with open(exon_list) as file:
        contents = file.readlines()
    all_exons = [line.strip().split('\t') for line in contents]
    with open(guide_list) as file:
        contents = file.readlines()
    all_guides = [line.strip().split('\t') for line in contents]
    columns = len(all_guides[0])
    for guide in all_guides[1:]:
        for exon in all_exons:
            if exon[0] == guide[1]:
                if guide[4] == '1':
                    cut_site = int(guide[3]) + 21
                elif guide[4] == '-1':
                    cut_site = int(guide[3]) + 9
                if cut_site >= int(exon[2]) and cut_site < int(exon[3])+1:
                    if len(guide) == columns:
                        guide.append(exon[6])
                        guide.append(exon[4])
                    #elif len(guide) == columns+1:
                    #    if guide[4] == '1':
                    #        if exon[6] > guide[columns]:
                    #            guide[columns] = exon[6]
                    #    if guide[4] == '-1':
                    #        if exon[6] > guide[columns]:
                    #            guide[columns] = exon[6]
    for guide in all_guides[1:]:
        if len(guide) == columns:
            guide.append('0')
            guide.append('0')
    all_guides[0].append('exon_score')
    all_guides[0].append('exon_dir')
    lines = []
    for i in all_guides:
        lines.append('\t'.join(i))
    with open(outfile, 'w') as file:
        file.write('\n'.join(lines))

def gc_content(guide_list, outfile):
    """
    add gc content (percent) column to guide list
    """
    with open(guide_list) as file:
        contents = file.readlines()
    all_guides = [line.strip().split('\t') for line in contents]
    for guide in all_guides[1:]:
        if guide[5] == '1':
            guide_seq = guide[8][4:24]
        elif guide[5] == '-1':
            guide_seq = guide[8][6:26]
        score = []
        for base in guide_seq:
            if base == 'G' or base == 'C':
                score.append(1)
        gc = sum(score) * 5
        guide.append(str(gc))
    all_guides[0].append('gc_percent')
    lines = []
    for i in all_guides:
        lines.append('\t'.join(i))
    with open(outfile, 'w') as file:
        file.write('\n'.join(lines))

if __name__ == '__main__':
    cli()