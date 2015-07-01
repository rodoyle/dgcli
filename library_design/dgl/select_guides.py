import os
import subprocess
import re
import math
import click

from dgl import library_constants
from operator import itemgetter

@click.group()
def cli():
    pass

@cli.command()
@click.argument('guide_list')
@click.argument('outfile')
@click.argument('n_guides')
@click.argument('spread', default='100')
def guide_selection(guide_list, outfile, n_guides, spread):
    with open(guide_list) as file:
        contents = file.readlines()
    all_guides = [line.strip().split('\t') for line in contents]
    #remove high off-targets
    left_guides = []
    for guide in all_guides[1:]:
        if float(guide[9]) == 0.0 and float(guide[10]) <=0.2:
            if float(guide[15]) <= 1.0:
                left_guides.append(guide)
    lines = ['\t'.join(all_guides[0])]
    for i in left_guides:
        lines.append('\t'.join(i))
    with open('approved_guides.txt', 'w') as file:
        file.write('\n'.join(lines))
    genes = set([guide[2] for guide in all_guides[1:]])
    black_swans = []
    for gene in genes:
        left_gene_guides = [guide for guide in left_guides if guide[2] == gene]
        #right_gene_guides = [guide for guide in right_guides if guide[1] == gene]
        cut_sites=[]
        bs = []
        doench_sorted = sorted(left_gene_guides, key=itemgetter(7))
        s = int(spread)
        for guide in reversed(doench_sorted):
            if len(cut_sites) < int(n_guides):
                if guide[5] == '1':
                    cut_site = int(guide[4]) + 21
                elif guide[5] == '-1':
                    cut_site = int(guide[4]) + 9
                if len(cut_sites) == 0:
                    bs.append(guide)
                    cut_sites.append([cut_site-s, cut_site+s])
                else:
                    collisions = 0
                    for cuts in cut_sites:
                        if cut_site >= cuts[0] and cut_site <= cuts[1]:
                            collisions += 1
                    if collisions == 0:
                        bs.append(guide)
                        cut_sites.append([cut_site-s, cut_site+s])
        black_swans = bs + black_swans
    lines = ['\t'.join(all_guides[0])]
    for i in black_swans:
        lines.append('\t'.join(i))
    with open(outfile, 'w') as file:
        file.write('\n'.join(lines))

if __name__ == '__main__':
    cli()


