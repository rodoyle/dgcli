import click
import subprocess
import os
import re
import yaml
import csv

from dgl import library_constants


@click.group()
def cli():
    pass

@cli.command()
@click.argument('fasta_in')
@click.argument('index_prefix')
@click.argument('mm')
@click.argument('chromosome_constant')
@click.argument('five_prime_trim', default=0)
def run_bowtie(fasta_in, index_prefix, mm, chromosome_constant, five_prime_trim):
    """
    use bowtie to locate
    potential offtarget scores for (nt)bp
    protospacer allowing for mm mismatches
    five_prime_trim allows truncation of nt from 5' end
    """
    chromosomes = getattr(library_constants, chromosome_constant)
    for chr in chromosomes:
        subprocess.call(
                ['bowtie', index_prefix + chr, # specify index
                 '-f', fasta_in + '.fa', # input is fasta
                 #'-a', # show all alignments
                 '-S', # output is sam
                 '-5', str(five_prime_trim), # bases to trim off 5' end
                 '-n', str(mm), # allow mismatches
                 '-v', str(mm), # allow mismatches
                 #'-p', '8', # parallel search threads
                 '-k', '1000', # allow 10 alignments per read
                 fasta_in + '_chr' + chr + '.sam'])

@cli.command()
@click.argument('biodata_prefix')
@click.argument('guide_list')
@click.argument('nt')
def parse_bowtie_sam(biodata_prefix, guide_list, nt):
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
    guide_dict = {guide[0]:guide[1] for guide in all_guides if len(guide) > 1}
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
            if offt_list.get(alignment[0]) == None:
                offt_list[alignment[0]] = []
            complement = {
                    'A':'T',
                    'T':'A',
                    'C':'G',
                    'G':'C'
                    }
            if strand[alignment[1]] == '-1':
                protospacer_seq = ''.join([complement[i] for i in reversed(alignment[9])])
            else:
                protospacer_seq = alignment[9]
            offt = {'guide':alignment[0],
                    'chr':chr,
                    'cut_site':int(alignment[3]),
                    'strand':strand[alignment[1]]}
            chr_pos = int(alignment[3])
            nt = len(alignment[9])
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
                    offt['pam'] = ''.join([complement[i] for i in reversed(
                        chr_fasta[chr_pos-4:chr_pos-1])])
                    offt['protospacer'] = ''.join([complement[i] for i in reversed(
                        chr_fasta[chr_pos-1:chr_pos+19])])
                    # find mismatch positions
                    mm = []
                    for nt in range(0,20,1):
                        if offt['protospacer'][nt] != guide_dict[alignment[0]][nt]:
                            mm.append(str(nt+1))
                    offt['mismatches'] = mm
                    offt_list[alignment[0]].append(offt)
    print yaml.dump({alignment[0]:offt_list})

@cli.command()
@click.argument('header')
def parse_guidebook(header):
    """
    parse guidebook output pasted into a header + _gb.txt file.
    output is a yaml-dump
    """
    with open(header + '_gb.txt') as file:
        contents = file.readlines()
    gbs = [line.strip().split('\t') for line in contents]
    gb_seqs = [gbs[line] for line in range(0,len(contents),3)]
    gb_pams = [gbs[line] for line in range(1,len(contents),3)]
    gb_infos = [gbs[line] for line in range(2,len(contents),3)]
    offt_list = []
    for n in range(0,len(gb_seqs),1):
        c1 = gb_infos[n][3].split('@ ')
        chr = c1[0].split('chr')[1]
        pos = int(c1[1].split('-')[0])
        offt = {'guide':header,
                'chr':chr,
                'cut_site':pos,
                'strand':'0',
                'protospacer':gb_seqs[n][0],
                'pam':gb_pams[n][0],
                'mismatches':gb_infos[n][1].split(',')}
        offt_list.append(offt)
    print yaml.dump({header:offt_list})

@cli.command()
@click.argument('header')
@click.argument('original_protospacer')
def parse_mit(header, original_protospacer):
    """
    parse mit csv output file.
    original protospacer = 20bp on-target sequence
    output is a yaml-dump
    """
    with open(header + '_offtargets.csv', 'rb') as csvfile:
        contents = csv.reader(csvfile, delimiter=',')
        mits = []
        for row in contents:
            mits.append(row)
    offt_list = []
    for mit in mits[1:]:
        mm=[]
        for nt in range(0,20,1):
            if original_protospacer[nt] != mit[4][nt+1]:
                mm.append(str(nt+1))
        offt = {'guide':header,
                'chr':mit[1].split('r')[1],
                'cut_site':int(mit[3])+15,
                'strand':mit[2],
                'protospacer':mit[4][0:20],
                'pam':mit[4][20:],
                'mismatches':mm}
        offt_list.append(offt)
    print yaml.dump({header:offt_list})

@cli.command()
def parse_guideseq():
    """
    parse guideseq.csv from xslx file.
    output is a yaml-dump
    """
    with open('guideseq.csv', 'rU') as csvfile:
        contents = csv.reader(csvfile, delimiter=',')
        guideseqs = []
        for row in contents:
            guideseqs.append(row)
    strand = {'+':'1',
              '-':'-1'}
    offt_list = {}
    for gs in guideseqs[1:]:
        if gs[7].find('site') > -1:
            sitename = gs[7].split('_site')
            gs_name = sitename[1] + '_' + sitename[0] + '_guideseq'
        else:
            gs_name = '1_' + gs[7] + '_guideseq'
        if offt_list.get(gs_name) == None:
            offt_list[gs_name] = []
        if gs[5] == '+':
            cut = 15
        elif gs[5] == '-':
            cut = 8
        mm=[]
        for nt in range(0,20,1):
            if gs[8][nt] != gs[9][nt]:
                mm.append(str(nt+1))
        offt = {'guide':gs_name,
                'chr':gs[0].split('r')[1],
                'cut_site':int(gs[1])+cut,
                'strand':strand[gs[5]],
                'protospacer':gs[9][0:20],
                'pam':gs[9][20:],
                'mismatches':mm,
                'guideseq_reads':gs[4]}
        offt_list[gs_name].append(offt)
    print yaml.dump(offt_list)

def compare_offt_outputs(guide_name, output_file1, output_file2):
    """
    compare two different off-target outputs
    in yaml format
    """
    with open(output_file1) as file:
        contents = file.read()
    offt_list = yaml.load(contents)

if __name__ == '__main__':
    cli()
