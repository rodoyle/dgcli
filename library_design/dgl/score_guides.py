import os
import subprocess
import re
import math
import click
import csv

from dgl import library_constants
from genome_reannotate import parse_gtf

COMPLEMENT = {'A':'T',
              'T':'A',
              'C':'G',
              'G':'C',
              'N':'N'
              }

@click.group()
def cli():
    pass

def doench(sequence):
    """
    Based on original Doench scoring algorithm from
    Doench et al 2014.
    """
    score = getattr(library_constants, 'INTERCEPT')
    # get the GC count
    gc_count = sequence[4:-6].count("G") + sequence[4:-6].count("C")
    # slice off the first 4 and last 3 and PAM to get protospacer ONLY
    gc_weight = getattr(library_constants, 'GC_HIGH')
    if gc_count <= 10:
        gc_weight = getattr(library_constants, 'GC_LOW')
    score += math.fabs(10 - gc_count) * gc_weight
    def get_weights(triplet):
        pos, model_bases, weight = triplet
        sub_sequence = sequence[pos:pos + len(model_bases)]
        if sub_sequence == model_bases:
            return weight
        return 0.0
    score += sum(map(get_weights, getattr(library_constants, 'DOENCH_DATA')))
    return 1.0 / (1.0 + math.exp(-score))

@cli.command()
@click.argument('guide_list')
@click.argument('outfile')
def on_target_score(guide_list, outfile):
    """
    calculate Doench score for guides and 
    append to guide list.
    also add header line to guide list
    """
    with open(guide_list) as file:
        contents = file.readlines()
    all_guides = [line.strip().split('\t') for line in contents]
    for guide in all_guides:
        if guide[5] == '-1':
            complement = {
            'A':'T',
            'T':'A',
            'C':'G',
            'G':'C'
            }
            guideseq = ''.join([complement[i] for i in reversed(guide[6])])
        else:
            guideseq = guide[6]
        score = doench(guideseq)
        guide.append("%.3f" % score)
    lines = []
    for i in all_guides:
        lines.append('\t'.join(i))
    with open(outfile, 'w') as file:
        file.write('guide_id\tgene_name\tgene_id\tchromosome\tabs_start\tstrand\tsequence\tdoench\n')
        file.write('\n'.join(lines))

def make_offtarget_fasta(guide_list, outfile, nt=15):
    """
    guide seq in guide_list[5]
    Reads in 30bp guides, and converts into
    (nt)bp guides fasta - default is 15 PAM-proximal for offtarget search
    """
    with open(guide_list) as file:
        contents = file.readlines()
    all_guides = [line.strip().split('\t') for line in contents]
    f = open(outfile, 'w')
    f.close()
    for guide in all_guides[1:]:
        if guide[5] == '-1':
            seq_end = int(nt)+6
            complement = {
                'A':'T',
                'T':'A',
                'C':'G',
                'G':'C'
                }
            guideseq = ''.join([complement[i] for i in reversed(guide[6][6:seq_end])])
        elif guide[5] == '1':
            seq_end = 24-int(nt)
            guideseq = guide[6][seq_end:24]
        with open(outfile, 'a') as file:
            file.write('>' + guide[0] + '\n')
            file.write(guideseq + '\n')

@cli.command()
@click.argument('fasta_prefix')
@click.argument('biodata')
@click.argument('chromosome_constant')
def build_bowtie_indexes(fasta_prefix, biodata, chromosome_constant):
    """
    build bowtie indexes for each
    chromosome
    """
    chromosomes = getattr(library_constants, chromosome_constant)
    for chr in chromosomes:
        fasta_file = biodata + fasta_prefix + chr + '.fasta'
        subprocess.call(
        ['bowtie-build',
         fasta_file,
         fasta_prefix + chr])

def offtarget_search_bowtie(fasta_in, sam_prefix, index_prefix, chromosome_constant, mm):
    """
    use bowtie to locate
    potential offtarget scores for (nt)bp (from fasta)
    protospacer allowing for (mm) mismatches
    """
    chromosomes = getattr(library_constants, chromosome_constant)
    for chr in chromosomes:
        subprocess.call(
                ['bowtie', index_prefix + chr, # specify index
                 '-f', fasta_in, # input is fasta
                 '-S', # output is sam
                 '-n', str(mm), # allow mismatches
                 '-p', '8', # parallel search threads
                 '-k', '10000', # allow 1000 alignments per read
                 '-y', # slow, comprehensive search
                 '-e', '1000', # combined phred score of mismatches to reject alignment
                 sam_prefix + chr + '.sam'])

def sam_extract_bowtie(biodata_prefix, nt):
    """
    find downstream 3bp from chromosome fasta
    remove non-PAM targets
    get mismatch positions from CIGAR string
    combine .sam files for each chromosome
    run script from sam file-containing directory
    """
    files = os.listdir('.')
    sam_files = [f for f in files if f.endswith('.sam')]
    offt_list=[]
    strand = {'0':'1',
              '256':'1',
              '16':'-1',
              '272':'-1'}
    COMPLEMENT = {'A':'T',
                  'T':'A',
                  'C':'G',
                  'G':'C'
                  }
    for sam in sam_files:
        chr1 = sam.split('.')
        chr2 = chr1[0].split('chr')
        chr = chr2[1]
        with open(biodata_prefix + chr + '.fasta') as file:
            sequence = file.read()
        chr_fasta = sequence.strip().split('\n')[1]
        with open(sam) as file:
            contents = file.readlines()
        all_alignments = [line.strip().split('\t') for line in contents if line.startswith('@')==False]
        alignments = [line for line in all_alignments if line[5]!='*']
        for alignment in alignments:
            offt = [alignment[0], chr, alignment[3], strand[alignment[1]]]
            chr_pos = int(alignment[3])
            nt = len(alignment[9])
            if alignment[1]=='0' or alignment[1]=='256':
                # if PAM site present downstream
                if re.findall("(?=([ACGT][AG]G))",
                              chr_fasta[chr_pos+nt-1:chr_pos+nt+2]) == [
                              chr_fasta[chr_pos+nt-1:chr_pos+nt+2]]:
                    # find mismatch positions
                    nt_mm = str()
                    mm_entry = alignment[12].split(':')
                    cigar = re.split(r"([A-Z])", mm_entry[2])
                    n = 0
                    for i in range(0,len(cigar)-1,2):
                        n = n + int(cigar[i]) + 1
                        offt.append(str(n+(20-int(nt))))
                    offt_list.append(offt)
            elif alignment[1]=='16' or alignment[1]=='272':
                # if PAM site present downstream (upstream)
                if re.findall("(?=(C[CT][ACGT]))",
                              chr_fasta[chr_pos-4:chr_pos-1]) == [
                              chr_fasta[chr_pos-4:chr_pos-1]]:
                    # find mismatch positions
                    nt_mm = str()
                    mm_entry = alignment[12].split(':')
                    cigar = re.split(r"([A-Z])", mm_entry[2])
                    n = 0
                    for i in range(0,len(cigar)-1,2):
                        n = n + int(cigar[i])+1
                        offt.append(str(21-n))
                    offt_list.append(offt)
    return offt_list

def parse_bowtie_sam(biodata, fasta_prefix, guide_list, nt):
    """
    find downstream 3bp from chromosome fasta (PAM)
    remove non-N[AG]G PAM-proximal targets
    get mismatch positions from CIGAR string
    combine .sam files for each chromosome
    run script from sam file-containing directory
    """
    files = os.listdir('.')
    sam_files = [f for f in files if f.endswith('.sam')]
    offt_list={}
    with open(guide_list) as file:
            contents = file.readlines()
    all_guides = [line.strip().split('\t') for line in contents]
    guide_dict = {}
    for guide in all_guides[1:]:
        if guide[5] == '1':
            guideseq = guide[6][4:24]
        elif guide[5] == '-1':
            guideseq = ''.join([COMPLEMENT[i] for i in reversed(guide[6][6:26])])
        guide_dict[guide[0]] = guideseq
    strand = {'0':'1',
              '256':'1',
              '16':'-1',
              '272':'-1'}
    for sam in sam_files:
        chr1 = sam.split('.')
        chr2 = chr1[0].split('chr')
        chr = chr2[1]
        with open(biodata + fasta_prefix + chr + '.fasta') as file:
            sequence = file.read()
        chr_fasta = sequence.strip().split('\n')[1]
        with open(sam) as file:
            contents = file.readlines()
        all_alignments = [line.strip().split('\t') for line in contents if line.startswith('@')==False]
        alignments = [line for line in all_alignments if line[5]!='*']
        for alignment in alignments:
            chr_pos = int(alignment[3])
            if offt_list.get(alignment[0]) == None:
                offt_list[alignment[0]] = []
            nt = len(alignment[9])
            if alignment[1]=='0' or alignment[1]=='256':
                cut_site = chr_pos+nt-6
            elif alignment[1]=='16' or alignment[1]=='272':
                cut_site = chr_pos+4
            offt = {'guide':alignment[0],
                    'chr':chr,
                    'cut_site':cut_site,
                    'strand':strand[alignment[1]]}
            if alignment[1]=='0' or alignment[1]=='256':
                # if PAM site present downstream
                if re.findall("(?=([ACGT][AG]G))",
                              chr_fasta[chr_pos+nt-1:chr_pos+nt+2]) == [
                              chr_fasta[chr_pos+nt-1:chr_pos+nt+2]]:
                    offt['pam'] = chr_fasta[chr_pos+nt-1:chr_pos+nt+2]
                    offt['protospacer'] = chr_fasta[chr_pos+nt-21:chr_pos+nt-1]
                    # find mismatch positions
                    mm = []
                    for nt in range(0,20,1):
                        if offt['protospacer'][nt] != guide_dict[alignment[0]][nt]:
                            mm.append(str(nt+1))
                    offt['mismatches'] = mm
                    offt_list[alignment[0]].append(offt)
            elif alignment[1]=='16' or alignment[1]=='272':
                # if PAM site present downstream (upstream)
                if re.findall("(?=(C[CT][ACGT]))",
                              chr_fasta[chr_pos-4:chr_pos-1]) == [
                              chr_fasta[chr_pos-4:chr_pos-1]]:
                    offt['pam'] = ''.join([COMPLEMENT[i] for i in reversed(
                        chr_fasta[chr_pos-4:chr_pos-1])])
                    offt['protospacer'] = ''.join([COMPLEMENT[i] for i in reversed(
                        chr_fasta[chr_pos-1:chr_pos+19])])
                    # find mismatch positions
                    mm = []
                    for nt in range(0,20,1):
                        if offt['protospacer'][nt] != guide_dict[alignment[0]][nt]:
                            mm.append(str(nt+1))
                    offt['mismatches'] = mm
                    offt_list[alignment[0]].append(offt)
    return offt_list

def mit_offtarget_score(offt_list):
    """
    Read in mismatch positions
    apply MIT scores to each off target hit
    """
    crispr_weights = getattr(library_constants, 'CRISPR_WEIGHTS')
    for guide in offt_list:
        for offt in offt_list[guide]:
            mismatch_count = len(offt['mismatches'])
            total_distance = 0.0
            if offt['pam'][1:3] == 'GG':
                score = 1.0
                score1 = 1.0
            elif offt['pam'][1:3] == 'AG':
                score = 0.2
                score1 = 0.2
            if mismatch_count == 0:
                offt['offt_score'] = score
            else:
                for n,mm in enumerate(offt['mismatches']):
                    if n > 0:
                        total_distance += (int(offt['mismatches'][n]) - int(offt['mismatches'][n - 1]))
                    score1 *= (1 - crispr_weights[int(mm)-1])
                mean_distance = 0.0 if mismatch_count <= 1 else (
                    total_distance / (mismatch_count - 1))
                score = score1 / ((19.0 - mean_distance) / 19.0 * 4.0 + 1.0)
                score /= mismatch_count * mismatch_count
                if mismatch_count == 1:
                    offt['offt_score'] = score1
                elif mismatch_count > 1:
                    offt['offt_score'] = score
    return offt_list

def is_coding_fast(biodata, embl_gtf, offt_list, nt=15):
    """
    uses Matt's colocator script to quickly seach fo overlapping regions
    """
    anns = parse_gtf(biodata + embl_gtf)
    cds_list=[]
    for cds in anns:
        if anns[cds][0] == 'CDS':
            if len(anns[cds][1]) < 3:
                cds_list.append(','.join([anns[cds][1],
                                         str(anns[cds][2]-1),
                                         str(anns[cds][3]),
                                         anns[cds][4],
                                         anns[cds][5]['gene_id']]))
    with open('cds_list.csv', 'w') as file:
        file.write(','.join(['chromosome_name',
                             'start',
                             'end',
                             'strand',
                             'label'
                             ]) + '\n')
        file.write('\n'.join(set(cds_list)))
    off_t_list = []
    offt_dict = {}
    for guide in offt_list:
        if offt_dict.get(guide) == None:
            offt_dict[guide] = {}
        for n,offt in enumerate(offt_list[guide]):
            offt_id = offt['guide']+'|'+str(n)
            offt_dict[guide][offt_id] = offt
            if offt['offt_score'] > 0.01:
                off_t_list.append(','.join([offt['chr'],
                                           str(offt['cut_site'] - 1),
                                           str(offt['cut_site'] + 1),
                                           offt['strand'],
                                           offt_id]))
    with open('offt_list.csv', 'w') as file:
        file.write(','.join(['chromosome_name',
                             'start',
                             'end',
                             'strand',
                             'label'
                             ]) + '\n')
        file.write('\n'.join(off_t_list))
    out_file = open('offt_colo.csv', 'w')
    subprocess.check_call(
        ['colocator',
         '--search', 'offt_list.csv',
         '--target', 'cds_list.csv',
         #'--output', 'offt_colo.csv',
         '--ignore-strand', '1', # only look at coordinates to judge overlap
         ],
         stdout=out_file)
    out_file.close()
    with open('offt_colo.csv', 'r') as file:
        contents = file.readlines()
    offtcs = [i.strip().split(',') for i in contents]
    for coding in offtcs:
        guide = coding[0].split('|')[0]
        offt_dict[guide][coding[0]]['coding'] = coding[1]
    hits = ['Off-Target_ID,Cut_Site,Strand,Protospacer,PAM,Mismatch_Positions,Off-Target_Score,Off-Target_Coding']
    for guide in offt_dict:
        for offt in offt_dict[guide]:
            if offt_dict[guide][offt].get('coding'):
                coding = offt_dict[guide][offt]['coding']
            else:
                coding = 'non-coding'
            mm='_'.join(offt_dict[guide][offt]['mismatches'])
            hits.append(','.join([offt,
                                  offt_dict[guide][offt]['chr'],
                                  str(offt_dict[guide][offt]['cut_site']),
                                  offt_dict[guide][offt]['strand'],
                                  offt_dict[guide][offt]['protospacer'],
                                  offt_dict[guide][offt]['pam'],
                                  mm,
                                  "%.6f" % offt_dict[guide][offt]['offt_score'],
                                  coding
                                  ]))
    with open('offt_hits.csv', 'w') as file:
        file.write('\n'.join(hits))
    return offt_dict

def count_offt_scores(guide_list, all_offt):
    with open(guide_list) as file:
        contents = file.readlines()
    all_guides = [line.strip().split('\t') for line in contents]
    for guide in all_guides[1:]:
        coding = 0
        genomic = 0
        self_off = 0
        gene_id = guide[0].split('_')[1]
        for offt in all_offt[guide[0]]:
            if all_offt[guide[0]][offt].get('coding'):
                coding += all_offt[guide[0]][offt]['offt_score']
                if all_offt[guide[0]][offt]['coding'] == gene_id:
                    self_off += all_offt[guide[0]][offt]['offt_score']
            else:
                genomic += all_offt[guide[0]][offt]['offt_score']
        guide.append("%.6f" % self_off)
        other_coding = coding - self_off
        guide.append("%.6f" % other_coding)
        guide.append("%.6f" % genomic)
    all_guides[0].append('self_offt')
    all_guides[0].append('coding_offt')
    all_guides[0].append('genomic_offt')
    return all_guides, all_offt

@cli.command()
@click.argument('guide_list')
@click.argument('outfile')
@click.argument('chromosome_constant')
@click.argument('fasta_prefix')
@click.argument('embl_gtf')
@click.argument('index_prefix')
@click.argument('biodata', default='/home/neil/biodata/')
@click.argument('nt', default='15')
def off_target_score(guide_list, outfile, chromosome_constant, fasta_prefix, embl_gtf, index_prefix, biodata, nt):
    """
    perform bowtie-based off-target scoring
    (build indexes using build_bowtie_indexes function)
    need to be in /Bowtie directory
    """
    make_offtarget_fasta(guide_list, '../guidelist.fasta')
    offtarget_search_bowtie('../guidelist.fasta', 'offt_alignment_chr',
                            biodata + index_prefix,
                            chromosome_constant, 
                            3)
    offt_list = parse_bowtie_sam(biodata, fasta_prefix, guide_list, nt)
    offt_scored = mit_offtarget_score(offt_list)
    del offt_list
    offt_coding = is_coding_fast(biodata, embl_gtf, offt_scored, nt)
    del offt_scored
    offt_guides = count_offt_scores(guide_list, offt_coding)
    lines = []
    for i in offt_guides[0]:
        lines.append('\t'.join(i))
    with open(outfile, 'w') as file:
        file.write('\n'.join(lines))

@cli.command()
@click.argument('guide_file')
@click.argument('index_prefix')
@click.argument('fasta_prefix')
@click.argument('embl_gtf')
@click.argument('biodata', default='/home/neil/biodata/')
def gecko_bowtie_search(guide_file, index_prefix, fasta_prefix, embl_gtf, biodata):
    # with open(guide_file) as file:
    #     contents = file.readlines()
    # guides = [i.strip().split(',') for i in contents]
    # with open('gecko_guides.fasta', 'w') as file:
    #     for guide in guides[1:]:
    #         file.write('>{}__{}\n'.format(guide[0], guide[1]))
    #         file.write(guide[2] + '\n')
    # with open('gecko_guidelist.txt', 'w') as file:
    #     file.write('0\t0\t0\t0\t0\t1\tNNNN0NNNNNN\n')
    #     for guide in guides[1:]:
    #        file.write('{}__{}\t0\t0\t0\t0\t1\tNNNN{}NNNNNN\n'.format(
    #            guide[0],
    #            guide[1],
    #            guide[2],
    #        ))
    # offtarget_search_bowtie('gecko_guides.fasta', # fasta in
    #                        'offt_alignment_chr', # sam prefix
    #                        biodata + index_prefix,
    #                        'LADY_CHROMOSOMES',
    #                        3)
    offt_list = parse_bowtie_sam(biodata, fasta_prefix, 'gecko_guidelist.txt', '20')
    offt_scored = mit_offtarget_score(offt_list)
    offt_coding = is_coding_fast(biodata, embl_gtf, offt_scored, nt=20)
    from pprint import pprint; import pdb; pdb.set_trace()
    # need to output self-scoring offtarget guides as csv

if __name__ == '__main__':
    cli()