import click
from random import randint


BASES = {
    0:'A',
    1:'C',
    2:'G',
    3:'T'
}

COMPLEMENT = {
    'A':'T',
    'T':'A',
    'C':'G',
    'G':'C'
}

def random_seq(length):
    sequence = str()
    for i in range(0,length-1,1):
        sequence = sequence + BASES[randint(0,3)]
    return sequence


@click.group()
def cli():
    pass

@cli.command()
@click.argument('guide_list')
#@click.argument('gene_list')
@click.argument('flanking_nt')
@click.argument('outfile_prefix')
def ready_guides(guide_list, flanking_nt, outfile_prefix):
    """
    Get 20bp forward protospacer sequence
    add additional required sequence
    flanking nt is a list: 5'Fx3'Fx5'Rx3'R
    output is 1x clone-able, and 1x information tab-delimited file
    """
    #with open(gene_list) as file:
    #    contents = file.read()
    #gene_lines = contents.strip().split('\r')
    #all_genes = [line.strip().split('\t') for line in gene_lines]
    #genes = {gene[1]:gene[0] for gene in all_genes}
    with open(guide_list) as file:
        contents = file.readlines()
    all_guides = [line.strip().split('\t') for line in contents]
    # Make Header Rows
    order_list = [['Guide_Number',
                   'Guide_ID',
                   '20bp_Protospacer',
                   'PAM',
                   'Oligo_Sequence']]
    info_list = [['Guide_Number',
                  'Guide_ID',
                  'Gene',
                  'Ensembl_ID',
                  'Chromosome',
                  'Position',
                  'Strand',
                  '20bp_Protospacer',
                  'PAM',
                  'On-Target_Doench',
                  'Off-Target_Self',
                  'Off-target_Coding',
                  'Off-target_Non-coding',
                  'Microhomology_Score',
                  'Out-of-Frame_Score',
                  'Exon_Score',
                  'Exon_Direction',
                  'CDS_Position',
                  'GC_percent',
                  'UUU_Occurences',
                  'Gs_after_PAM,',
                  'Gs_at_Start']]
    flanking = flanking_nt.split('x')
    for n,guide in enumerate(all_guides[1:]):
        if guide[5] == '-1':
            guideseq = ''.join([COMPLEMENT[i] for i in reversed(guide[6][6:26])])
            pam = ''.join([COMPLEMENT[i] for i in reversed(guide[6][3:6])])
            f_order = flanking[0] + \
                      guideseq + \
                      flanking[1]
            r_order = flanking[2] + \
                      ''.join([COMPLEMENT[i] for i in reversed(guideseq)]) + \
                      flanking[3]
        else:
            guideseq = guide[6][4:24]
            pam = guide[6][24:27]
            f_order = flanking[0] + \
                      guideseq + \
                      flanking[1]
            r_order = flanking[2] + \
                      ''.join([COMPLEMENT[i] for i in reversed(guideseq)]) + \
                      flanking[3]
        order_list.append([str(n+1),
                           guide[0],
                           guideseq,
                           pam,
                           f_order])
        order_list.append([str(n+1),
                           guide[0],
                           guideseq,
                           pam,
                           r_order])
        info_list.append([str(n+1)] + \
                         guide[0:6] + \
                         [guideseq, pam] + \
                         guide[7:15] + \
                         ["%.3f" % float(guide[15])] + \
                         guide[16:])
    lines = []
    for i in order_list:
        lines.append('\t'.join(i))
    with open(outfile_prefix + '_oligo_order.txt', 'w') as file:
        file.write('\n'.join(lines))
    lines = []
    for i in info_list:
        lines.append('\t'.join(i))
    with open(outfile_prefix + '_guide_info.txt', 'w') as file:
        file.write('\n'.join(lines))

if __name__ == '__main__':
    cli()