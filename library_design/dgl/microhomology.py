from math import exp
import click

from dgl import library_constants

@click.group()
def cli():
    pass

def microhomology_score(sequence):
    """
    runs script from Bae et al. 2014
    which takes a cut-site flanking sequence of 80bp and cut site location
    and returns [microhomology score, out-of-frame score, list_output]
    """
    seq = sequence
    length_weight=20.0
    left=30   # Insert the position of the cut site.
    right=len(seq)-left
    dups = []
    mh_results = []
    # nt positions from 39-2
    for k in range(2, left)[::-1]:
        # nt positions from 40-(42-79)
        for j in range(left, left+right-k+1):
            # range 0-(2-39)
            for i in range(0, left-k+1):
                # if sequences are homologous
                if seq[i:i+k] == seq[j:j+k]:
                    length=j-i
                    dups.append([seq[i:i+k], # homologous sequence
                                 i, # left 5' edge
                                 i+k, # left 3' edge
                                 j, # right 5' edge
                                 j+k, # right 3' edge
                                 length # left 5' to right 5'
                                 ])
    if len(dups) > 0:
        sum_score_3 = 0
        sum_score_not_3 = 0
        for i in range(len(dups)):
            n = 0
            score_3 = 0
            score_not_3 = 0
            scrap = dups[i][0]
            left_start = dups[i][1]
            left_end = dups[i][2]
            right_start = dups[i][3]
            right_end =dups[i][4]
            length = dups[i][5]
            # range 0-i
            for j in range(i):
                left_start_ref = dups[j][1]
                left_end_ref = dups[j][2]
                right_start_ref = dups[j][3]
                right_end_ref = dups[j][4]
                if i > 0:
                    import pdb; pdb.set_trace
                    # checks for duplicates
                    # if left_ref contains left and right_ref contains right
                    if (left_start >= left_start_ref) and \
                       (left_end <= left_end_ref) and \
                       (right_start >= right_start_ref) and \
                       (right_end <= right_end_ref):
                        # if left and right are the same length relative to ref
                        if (left_start - left_start_ref) == \
                           (right_start - right_start_ref) and \
                           (left_end - left_end_ref) == \
                           (right_end - right_end_ref):
                            n += 1
                        else:
                            pass
            if n == 0:
                length_factor = round(1/exp(length/length_weight),3)
                num_GC = scrap.count('G') + scrap.count('C')
                # if in-frame deletion
                if (length % 3) == 0:
                    score_3 = 100 * length_factor * ((len(scrap)-num_GC) + (num_GC*2))

                elif (length % 3) != 0:
                    score_not_3 = 100 * length_factor * ((len(scrap)-num_GC) + (num_GC*2))
                mh_results.append([seq[0:left_end] +\
                                  '-'*length + \
                                  seq[right_end:],
                                  scrap,
                                  str(length),
                                  str(100 * length_factor * ((len(scrap)-num_GC)+(num_GC*2)))
                                  ])
                sum_score_3 += score_3
                sum_score_not_3 += score_not_3
    mh_score = sum_score_3+sum_score_not_3
    oof_score = (sum_score_not_3)*100/(sum_score_3+sum_score_not_3)
    return ["%.1f" % mh_score, "%.2f" % oof_score, mh_results]

@cli.command()
@click.argument('guide_list')
@click.argument('biodata_prefix')
@click.argument('outfile')
def microhomology(guide_list, biodata_prefix, outfile):
    """
    reads in guidelist and adds columns for microhomology score and
    out-of-frame score. Uses 80bp cut-site-flanking sequence
    """
    with open(guide_list) as file:
        contents = file.readlines()
    all_guides = [line.strip().split('\t') for line in contents]
    chromosomes = set([g[3] for g in all_guides[1:]])
    mh_scored = []
    for chr in chromosomes:
        chr_guides = [g for g in all_guides if g[3] == chr]
        with open(biodata_prefix + chr + '.fasta') as file:
            sequence = file.read()
        chr_fasta = sequence.strip().split('\n')[1]
        for guide in chr_guides:
            if guide[5] == '1':
                cut_site = int(guide[4]) + 21
            elif guide[5] == '-1':
                cut_site = int(guide[4]) + 9
            scores = microhomology_score(chr_fasta[cut_site-41:cut_site+39])
            guide = guide + [scores[1]]
            mh_scored.append(guide)
    all_guides[0].append('out-of-frame_score')
    lines = ['\t'.join(all_guides[0])]
    for i in mh_scored:
        lines.append('\t'.join(i))
    with open(outfile, 'w') as file:
        file.write('\n'.join(lines))

if __name__ == '__main__':
    cli()