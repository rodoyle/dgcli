import subprocess
import os
import click
import yaml

from Bio import SeqIO
from dgl import library_constants

@click.group()
def cli():
    pass

@cli.command()
@click.argument('embl_gtf')
def parse_gtf(embl_gtf):
    """
    organise embl gtf into genes
    output in yaml format
    """
    with open(embl_gtf) as file:
        contents = file.readlines()
    anns = [line.strip().split('\t') for line in contents if line.startswith('#') == False]
    strand = {'+':'1',
              '-':'-1'}
    gene_anns = [i for i in anns if i[2] == 'gene']
    gene_dict = {}
    for gene in gene_anns:
        info_parsed = gene[8].split('; ')
        info_dict = {}
        for info in info_parsed:
            info_dict[info.split(' "')[0]] = info.split('"')[1]
        gene_dict[info_dict['gene_id']] = [gene[2],
                                           gene[0],
                                           int(gene[3]),
                                           int(gene[4]),
                                           strand[gene[6]],
                                           info_dict]
    for ann_type in ['transcript', 'CDS', 'exon']:
        gene_anns = [i for i in anns if i[2] == ann_type]
        for anno in gene_anns:
            info_parsed = anno[8].split('; ')
            info_dict = {}
            for info in info_parsed:
                info_dict[info.split(' "')[0]] = info.split('"')[1]
            gene_dict[info_dict['gene_id']].append([anno[2],
                                                    anno[0],
                                                    int(anno[3]),
                                                    int(anno[4]),
                                                    strand[anno[6]],
                                                    info_dict])
    import pdb; pdb.set_trace()
    with open('Rattus_norvegicus.Rnor_5.0.79.yaml', 'w') as file:
        file.write(yaml.dump(gene_dict))
    #print yaml.dump(gene_dict)



@cli.command()
@click.argument('embl_gtf')
@click.argument('fasta_prefix')
@click.argument('ann_type')
@click.argument('nt')
@click.argument('biodata')
@click.argument('chromosome_constant')
def get_ann_seqs_gtf(embl_gtf, fasta_prefix, ann_type, nt, biodata, chromosome_constant):
    """
    Use to prepare sequences for Subread aligner;
    head_tail argument makes head and tail files of
    first and last nt bp of each feature to map.
    Tail is reverse complement of end of annotation to get
    the end position
    """
    chromosomes = getattr(library_constants, chromosome_constant)
    with open(embl_gtf) as file:
        contents = file.readlines()
    anns = [line.strip().split('\t') for line in contents]
    f=open('too_short_' + ann_type + '.txt', 'w')
    f.close()
    complement = {
            'A':'T',
            'T':'A',
            'C':'G',
            'G':'C',
            'a':'T',
            't':'A',
            'c':'G',
            'g':'C',
            'N':'N'
            }
    for chr in chromosomes:
        f=open('head_' + fasta_prefix + chr + '.fasta', 'w')
        f.close()
        f=open('tail_' + fasta_prefix + chr + '.fasta', 'w')
        f.close()
        with open(biodata + chr + '.fasta') as file:
            sequence = file.read()
        chr_fasta = sequence.strip().split('\n')[1]
        chr_anns = [line.strip().split('\t') for line in contents if line.startswith(chr + '\t')]
        for n,genome_feature in enumerate(chr_anns):
            if genome_feature[2] == ann_type:
                if int(genome_feature[4]) - int(genome_feature[3]) +1 < 20:
                    with open('too_short_' + ann_type + '.txt', 'a') as file:
                        file.write('\t'.join(genome_feature) + '\n')
                else:
                    prefix = '>' + str(n) + '|'
                    if int(genome_feature[4]) - int(genome_feature[3]) +1 < int(nt):
                        with open('head_' + fasta_prefix + chr + '.fasta', 'a') as file:
                            file.write(prefix + '|'.join(genome_feature) + '\n')
                            file.write(chr_fasta[int(genome_feature[3])-1:int(genome_feature[4])] + '\n')
                        with open('tail_' + fasta_prefix + chr + '.fasta', 'a') as file:
                            file.write(prefix + '|'.join(genome_feature) + '\n')
                            file.write(''.join([
                                       complement[i] for i in reversed(
                                        chr_fasta[int(genome_feature[3])-1:int(genome_feature[4])])]) + '\n')
                    else:
                        with open('head_' + fasta_prefix + chr + '.fasta', 'a') as file:
                            file.write(prefix + '|'.join(genome_feature) + '\n')
                            file.write(chr_fasta[int(genome_feature[3])-1:int(genome_feature[3])+int(nt)-1] + '\n')
                        with open('tail_' + fasta_prefix + chr + '.fasta', 'a') as file:
                            file.write(prefix + '|'.join(genome_feature) + '\n')
                            file.write(''.join([
                                       complement[i] for i in reversed(
                                        chr_fasta[int(genome_feature[4])-int(nt):int(genome_feature[4])])]) + '\n')

@cli.command()
@click.argument('fasta_prefix')
@click.argument('biodata')
@click.argument('chromosomes')
def build_subread_indexes(fasta_prefix, biodata, chromosomes):
    """
    build subread aligner indexes for each
    chromosome
    """
    for chr in chromosomes:
        fasta_file = biodata + fasta_prefix + chr + '.fasta'
        subprocess.call(
        ['subread-buildindex',
         fasta_file,
         '-o', fasta_prefix + chr])



@cli.command()
@click.argument('fasta_prefix')
@click.argument('index_prefix')
@click.argument('chromosomes')
def run_subread_aligner(fasta_prefix, index_prefix, chromosomes):
    """
    align feature heads or tails to indexed genome using
    Subread aligner
    """
    for chr in chromosomes:
        subprocess.call(
            ['subread-align',
             '-T', '5',
             '-B', '1',
             '-i', index_prefix + chr,
             '-r', 'head_' + fasta_prefix + chr + '.fasta',
             '-o', 'head_' + fasta_prefix + chr + '.sam'])
        subprocess.call(
            ['subread-align',
             '-T', '5',
             '-B', '1',
             '-i', index_prefix + chr,
             '-r', 'tail_' + fasta_prefix + chr + '.fasta',
             '-o', 'tail_' + fasta_prefix + chr + '.sam'])

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

@cli.command()
@click.argument('fasta_in')
@click.argument('index_prefix')
@click.argument('mm')
@click.argument('chromosome_constant')
@click.argument('five_prime_trim', default='0')
def run_bowtie(fasta_in, index_prefix, mm, chromosome_constant, five_prime_trim):
    """
    use bowtie to locate
    potential offtarget scores for (nt)bp
    protospacer allowing for mm=3 mismatches
    """
    chromosomes = getattr(library_constants, chromosome_constant)
    for chr in chromosomes:
        subprocess.call(
                ['bowtie', index_prefix + chr, # specify index
                 '-f', fasta_in + chr + '.fasta', # input is fasta
                 '-S', # output is sam
                 '-n', str(mm), # allow mismatches
                 '-p', '8', # parallel search threads
                 '-k', '10', # allow 10 alignments per read
                 fasta_in + chr + '.sam'])

@cli.command()
@click.argument('file_name')
@click.argument('outfile')
@click.argument('chromosome_constant')
def confirm_features_bowtie(file_name, outfile, chromosome_constant):
    """
    file_name = 'fasta_prefix'
    confirm mapped coordinates using bowtie output 
    by checking if feature length is within 20% of expected
    length. Output as 3 lists, one containing newly-mapped
    annotations, one for non-mapped at start or end and one
    for incorrect lenth mapping
    """
    new_anns = []
    f = open(outfile, 'w')
    f.close
    wrong_length = []
    f = open('wrong_length.txt', 'w')
    f.close
    not_mapped = []
    f = open('not_mapped.txt', 'w')
    f.close
    multi_aligned = []
    f = open('multi_aligned.txt', 'w')
    f.close
    chromosomes = getattr(library_constants, chromosome_constant)
    STRAND = {'+':'1',
              '-':'-1'}
    for chr in chromosomes:
        with open('head_' + file_name + chr + '.sam', 'r') as file:
            contents = file.readlines()
        sam_lines = [line.strip().split('\t') for line in contents \
                     if line.startswith("@") == False]
        head_dict = {}
        for i in sam_lines:
            if i[1]=='0':
                if head_dict.get(i[0]):
                    head_dict[i[0]].append('multi')
                else:
                    head_dict[i[0]] = [i[3], i[4]]
        with open('tail_' + file_name + chr + '.sam', 'r') as file:
            contents = file.readlines()
        sam_lines = [line.strip().split('\t') for line in contents \
                     if line.startswith("@") == False]
        tail_dict = {}
        for i in sam_lines:
            if i[1]=='16':
                if tail_dict.get(i[0]):
                    tail_dict[i[0]].append('multi')
                else:
                    tail_dict[i[0]] = [i[3],len(i[9]), i[4]]
        for feature_id in head_dict:
            if len(head_dict[feature_id]) > 2:
                if head_dict[feature_id][2] == 'multi':
                    multi_aligned.append(['head',feature_id])
            else:
                id_parsed = feature_id.split('|')
                if tail_dict.get(feature_id):
                    if len(head_dict[feature_id]) > 3:
                        if tail_dict[feature_id][3] == 'multi':
                            multi_aligned.append(['tail',feature_id])
                    else:
                        length = (int(tail_dict[feature_id][0])+tail_dict[feature_id][1])-int(head_dict[feature_id][0])
                        exp_length = int(id_parsed[5])-int(id_parsed[4])
                        if float(length)*0.8 < float(exp_length) and float(length)*1.2 > float(exp_length):
                            new_anns.append([feature_id, 
                                             head_dict[feature_id][0], 
                                             str(int(tail_dict[feature_id][0]) + int(tail_dict[feature_id][1]) -1), 
                                             STRAND[id_parsed[7]]])
                        else:
                            difference = float(length)/float(exp_length)
                            wrong_length.append([feature_id, 
                                                 head_dict[feature_id][0], 
                                                 tail_dict[feature_id][0], 
                                                 "%.2f" % difference])
                else:
                    not_mapped.append(['tail',feature_id])
        for feature_id in tail_dict:
            if head_dict.get(feature_id):
                continue
            else:
                not_mapped.append(['head',feature_id])
    lines = []
    for i in new_anns:
        lines.append('\t'.join(i))
    with open(outfile, 'a') as file:
        file.write('\n'.join(lines))
        file.write('\n')
    lines = []
    for i in not_mapped:
        lines.append('\t'.join(i))
    with open('not_mapped.txt', 'a') as file:
        file.write('\n'.join(lines))
        file.write('\n')
    lines = []
    for i in wrong_length:
        lines.append('\t'.join(i))
    with open('wrong_length.txt', 'a') as file:
        file.write('\n'.join(lines))
        file.write('\n')
    lines = []
    for i in multi_aligned:
        lines.append('\t'.join(i))
    with open('multi_aligned.txt', 'a') as file:
        file.write('\n'.join(lines))
        file.write('\n')

@cli.command()
@click.argument('fasta_prefix')
@click.argument('biodata')
@click.argument('chromosome_constant')
def build_bwa_mem_indexes(fasta_prefix, biodata, chromosome_constant):
    """
    build indexes for bwa-mem
    """
    chromosomes = getattr(library_constants, chromosome_constant)
    for chr in chromosomes:
        subprocess.call([
            'bwa', 'index',
            '-p', fasta_prefix + chr, # prefix of output
            '-a', 'is', # type of algorithm
            biodata + fasta_prefix + chr + '.fasta' # fasta index file
            ])

@cli.command()
@click.argument('index_prefix')
@click.argument('fasta_prefix')
@click.argument('chromosome_constant')
@click.argument('fasta_number', default='many')
def run_bwa_mem(index_prefix, fasta_prefix, chromosome_constant, fasta_number):
    """
    run bwa-mem
    """
    chromosomes = getattr(library_constants, chromosome_constant)
    if fasta_number != '1':
        for chr in chromosomes:
            head_outfile = open('head_' + fasta_prefix + chr + '.sam', 'w')
            tail_outfile = open('tail_' + fasta_prefix + chr + '.sam', 'w')
            subprocess.check_call([
                'bwa', 'mem',
                '-k', '10', # minimum seed length
                index_prefix + chr,
                'head_' + fasta_prefix + chr + '.fasta'], 
                stdout=head_outfile)
            subprocess.check_call([
                'bwa', 'mem',
                '-k', '10', # minimum seed length
                index_prefix + chr,
                'tail_' + fasta_prefix + chr + '.fasta'], 
                stdout=tail_outfile)
    else:
        for chr in chromosomes:
            head_outfile = open('head_' + fasta_prefix + '.sam', 'w')
            tail_outfile = open('tail_' + fasta_prefix + '.sam', 'w')
            subprocess.check_call([
                'bwa', 'mem',
                '-k', '10', # minimum seed length
                index_prefix + chr,
                'head_' + fasta_prefix + '.fasta'], 
                stdout=head_outfile)
            subprocess.check_call([
                'bwa', 'mem',
                '-k', '10', # minimum seed length
                index_prefix + chr,
                'tail_' + fasta_prefix + '.fasta'], 
                stdout=tail_outfile)

@cli.command()
@click.argument('file_in')
@click.argument('fasta_prefix')
@click.argument('nt')
@click.argument('biodata')
@click.argument('chromosome_constant')
@click.argument('id_column', default=0)
def recycle_fasta(file_in, fasta_prefix, nt, biodata, chromosome_constant, id_column):
    """
    get fasta sequences for next round
    """
    chromosomes = getattr(library_constants, chromosome_constant)
    with open(file_in, 'r') as file:
        contents = file.readlines()
    useless = [line.strip().split('\t') for line in contents]
    if id_column == 0:
        identifiers = [i.split('|') for i in useless]
    else:
        identifiers = [i[int(id_column)].split('|') for i in useless]
    complement = {
            'A':'T',
            'T':'A',
            'C':'G',
            'G':'C',
            'a':'T',
            't':'A',
            'c':'G',
            'g':'C',
            'N':'N'
            }
    for chr in chromosomes:
        f=open('head_' + fasta_prefix + chr + '.fasta', 'w')
        f.close()
        f=open('tail_' + fasta_prefix + chr + '.fasta', 'w')
        f.close()
        with open(biodata + chr + '.fasta') as file:
            sequence = file.read()
        chr_fasta = sequence.strip().split('\n')[1]
        chr_anns = [i for i in identifiers if i[1]==chr]
        for genome_feature in chr_anns:
            if int(genome_feature[5]) - int(genome_feature[4]) +1 < int(nt):
                with open('head_' + fasta_prefix + chr + '.fasta', 'a') as file:
                    file.write('>' + '|'.join(genome_feature) + '\n')
                    file.write(chr_fasta[int(genome_feature[4])-1:int(genome_feature[5])] + '\n')
                with open('tail_' + fasta_prefix + chr + '.fasta', 'a') as file:
                    file.write('>' + '|'.join(genome_feature) + '\n')
                    file.write(''.join([
                               complement[i] for i in reversed(
                                chr_fasta[int(genome_feature[4])-1:int(genome_feature[5])])]) + '\n')
            else:
                with open('head_' + fasta_prefix + chr + '.fasta', 'a') as file:
                    file.write('>' + '|'.join(genome_feature) + '\n')
                    file.write(chr_fasta[int(genome_feature[4])-1:int(genome_feature[4])+int(nt)-1] + '\n')
                with open('tail_' + fasta_prefix + chr + '.fasta', 'a') as file:
                    file.write('>' + '|'.join(genome_feature) + '\n')
                    file.write(''.join([
                               complement[i] for i in reversed(
                                chr_fasta[int(genome_feature[5])-int(nt):int(genome_feature[5])])]) + '\n')

@cli.command()
@click.argument('list_input')
@click.argument('original_fasta')
@click.argument('outfile')
def get_annotations_back(list_input, original_fasta, outfile):
    with open(list_input, 'r') as file:
        contents = file.readlines()
    useless = [line.strip().split('\t') for line in contents]
    if len(useless[0]) == 2:
        n=1 # column with header data
    else:
        n=0
    with open(original_fasta, 'r') as file:
        contents = file.read()
    fastas = contents.strip().split('\n')
    headers = []
    for i in range(0,len(fastas),2):
        headers.append(fastas[i].split('>')[1])
    multi_anns = []
    for head in headers:
        for multi in useless:
            if head.startswith(multi[n]):
                multi_anns.append(head.split('|')[1:])
    lines = []
    for i in multi_anns:
        lines.append('\t'.join(i))
    with open(outfile, 'a') as file:
        file.write('\n'.join(lines))




if __name__ == '__main__':
    cli()