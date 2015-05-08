import os
import subprocess
import re
import math
import click

from dgl import library_constants

@click.group()
def cli():
    pass

def doench(sequence):
    """
    Based on original Doench scoring algorithm from
    Doench et al 2014.
    """
    score = getattr(library_constants, INTERCEPT)
    # get the GC count
    gc_count = sequence[4:-6].count("G") + sequence[4:-6].count("C")
    # slice off the first 4 and last 3 and PAM to get protospacer ONLY
    gc_weight = getattr(library_constants, GC_HIGH)
    if gc_count <= 10:
        gc_weight = getattr(library_constants, GC_LOW)
    score += math.fabs(10 - gc_count) * gc_weight
    def get_weights(triplet):
        pos, model_bases, weight = triplet
        sub_sequence = sequence[pos:pos + len(model_bases)]
        if sub_sequence == model_bases:
            return weight
        return 0.0
    score += sum(map(get_weights, getattr(library_constants, DOENCH_DATA)))
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
        if guide[4] == '-1':
            complement = {
            'A':'T',
            'T':'A',
            'C':'G',
            'G':'C'
            }
            guideseq = ''.join([complement[i] for i in reversed(guide[5])])
        else:
            guideseq = guide[5]
        score = doench(guideseq)
        guide.append("%.3f" % score)
    lines = []
    for i in all_guides:
        lines.append('\t'.join(i))
    with open(outfile, 'w') as file:
        file.write('guide_id\tgene\tchromosome\tabs_start\tstrand\tsequence\tdoench\n')
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
        if guide[4] == '-1':
            seq_end = int(nt)+6
            complement = {
                'A':'T',
                'T':'A',
                'C':'G',
                'G':'C'
                }
            guideseq = ''.join([complement[i] for i in reversed(guide[5][6:seq_end])])
        elif guide[4] == '1':
            seq_end = 24-int(nt)
            guideseq = guide[5][seq_end:24]
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
                 '-k', '10', # allow 10 alignments per read
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

def mit_offtarget_score(offt_list):
    """
    Read in mismatch positions
    apply MIT scores to each off target hit
    """
    all_offt = offt_list
    offt_scored = []
    for offt in all_offt:
        score = 1.0
        if len(offt) == 4:
            offt_scored.append(offt[0:4] + ['1.000'])
        for mm in offt[4:]:
            score *= 1-CRISPR_WEIGHTS[int(mm)-1]
            if mm == offt[len(offt)-1]:
                offt_scored.append(offt[0:4] + ["%.3f" % score])
    return offt_scored

def is_coding_chr(cds_list, offt_list, chromosome_constant, nt=15):
    """
    determine if off-target hits cut in CDS regions
    appends each chromosome to offt outfile as this takes a while...
    """
    with open(cds_list) as file:
        contents = file.readlines()
    all = [line.strip().split('\t') for line in contents if line.startswith('#')==False]
    all_cdss = [anno for anno in all if anno[2] == 'CDS']
    del all
    all_offt = offt_list
    columns = len(all_offt[0])
    f = open('off_target_hitlist.txt', 'w')
    f.close()
    offt_coding = []
    chromosomes = getattr(library_constants, chromosome_constant)
    for chr in chromosomes:
        chr_cdss = [cds for cds in all_cdss if cds[0] == chr]
        chr_offt = [offt for offt in all_offt if offt[1] == chr]
        for offt in chr_offt:
            if offt[3] == '1':
                cut_site = int(offt[2]) + int(nt)-5
            elif offt[3] == '-1':
                cut_site = int(offt[2]) + 5
            for cds in chr_cdss:
                if len(offt) > columns:
                    break
                elif len(offt) == columns:
                    if cut_site >= int(cds[3]):
                        if cut_site <= int(cds[4]):
                            offt.append('1')
        for offt in chr_offt:
            if len(offt) == columns:
                offt.append('0')
            offt_coding.append(offt)
        lines = []
        for i in chr_offt:
            lines.append('\t'.join(i))
        with open('off_target_hitlist.txt', 'a') as file:
            file.write('\n'.join(lines))
            file.write('\n')
    return offt_coding

def count_offt_scores(guide_list, all_offt):
    with open(guide_list) as file:
        contents = file.readlines()
    all_guides = [line.strip().split('\t') for line in contents]
    for guide in all_guides[1:]:
        coding = 0
        genomic = 0
        guide_offt = [offt for offt in all_offt if offt[0] == guide[0]]
        for offt in guide_offt:
            if float(offt[4]) > 0.8 and offt[5] == '1':
                coding += 1
            elif float(offt[4]) > 0.8 and offt[5] == '0':
                genomic += 1
        guide.append(str(coding))
        guide.append(str(genomic))
    all_guides[0].append('coding_offt')
    all_guides[0].append('genomic_offt')
    return all_guides

@cli.command()
@click.argument('guide_list')
@click.argument('outfile')
@click.argument('chromosome_constant')
@click.argument('biodata', default='/home/neil/biodata/')
@click.argument('nt', default='15')
def off_target_score(guide_list, outfile, chromosome_constant, biodata, nt):
    """
    perform bowtie-based off-target scoring
    (build indexes using build_bowtie_indexes function)
    need to be in /Bowtie directory
    """
    make_offtarget_fasta(guide_list, '../guidelist.fasta')
    offtarget_search_bowtie('../guidelist.fasta', 'offt_alignment_chr',
                            biodata + 'bowtie_index/mouse_GRCm38/GRCm38_chr',
                            chromosome_constant, 
                            2)
    offt_list = sam_extract_bowtie(biodata + 'genomes_fasta/mouse_GRCm38/GRCm38_chr',nt)
    offt_scored = mit_offtarget_score(offt_list)
    del offt_list
    offt_coding = is_coding_chr(biodata + 'embl_annotations/Mus_musculus.GRCm38.79.gtf',
                                offt_scored,
                                chromosome_constant)
    del offt_scored
    offt_guides = count_offt_scores(guide_list, offt_coding)
    lines = []
    for i in offt_guides:
        lines.append('\t'.join(i))
    with open(outfile, 'w') as file:
        file.write('\n'.join(lines))

if __name__ == '__main__':
    cli()