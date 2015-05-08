import subprocess
import os
import click
import yaml

from dgl import library_constants
from Bio import AlignIO

@click.group()
def cli():
    pass

#@cli.command()
#@click.argument('embl_gtf')
def parse_gtf(embl_gtf):
    """
    organise embl gtf into genes
    output in yaml format
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
            gene_dict[info_dict['gene_id'] + '|' +
                      info_dict['transcript_id']] = ann_info
        elif ann[2] == 'exon':
            gene_dict[info_dict['gene_id'] + '|' +
                      info_dict['exon_id']] = ann_info
        elif ann[2] == 'CDS':
            gene_dict[info_dict['gene_id'] + '|' +
                      info_dict['protein_id'] + '|' +
                      info_dict['exon_number']] = ann_info
    return gene_dict

@cli.command()
@click.argument('embl_gtf')
@click.argument('fasta_prefix')
@click.argument('nt')
@click.argument('type')
@click.argument('biodata')
@click.argument('chromosome_constant')
def get_ann_seqs_fasta(embl_gtf, fasta_prefix, nt, type, biodata, chromosome_constant):
    """
    get annotation sequences and output as fasta
    nt = range of nts (a-b)
    type = gene, transcript, exon or CDS
    """
    chromosomes = getattr(library_constants, chromosome_constant)
    if embl_gtf.endswith('.gtf'):
        ann_dict = parse_gtf(embl_gtf)
    elif embl_gtf.endswith('.yaml'):
        with open(embl_gtf, 'r') as file:
            contents = file.read()
        ann_dict = yaml.load(contents)
    for chr in chromosomes:
        f=open(fasta_prefix + chr + '.fasta', 'w')
        f.close()
        with open(biodata + chr + '.fasta') as file:
            sequence = file.read()
        chr_fasta = sequence.strip().split('\n')[1]
        for gene_id in ann_dict:
            if ann_dict[gene_id][1] == chr:
                if ann_dict[gene_id][0] == type:
                    length = ann_dict[gene_id][3] - ann_dict[gene_id][2]
                    nt_range = nt.split('-')
                    if int(nt_range[0]) <= length:
                        if int(nt_range[1]) >= length:
                            with open(fasta_prefix + chr + '.fasta', 'a') as file:
                                file.write('>' + gene_id + '\n')
                                file.write(
                                    chr_fasta[ann_dict[gene_id][2]-1:
                                              ann_dict[gene_id][3]] + '\n')

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
            outfile = open(fasta_prefix + chr + '.sam', 'w')
            subprocess.check_call([
                'bwa', 'mem',
                '-k', '10', # minimum seed length
                '-a', # show all alignments
                index_prefix + chr,
                fasta_prefix + chr + '.fasta'], 
                stdout=outfile)
    else:
        for chr in chromosomes:
            outfile = open(fasta_prefix + chr + '.sam', 'w')
            subprocess.check_call([
                'bwa', 'mem',
                '-k', '10', # minimum seed length
                '-a', # show all alignments
                index_prefix + chr,
                fasta_prefix + '.fasta'], 
                stdout=outfile)

@cli.command()
@click.argument('fasta_in')
@click.argument('index_prefix')
@click.argument('mm')
@click.argument('chromosome_constant')
@click.argument('fasta_number')
@click.argument('five_prime_trim', default='0')
def run_bowtie(fasta_in, index_prefix, mm, chromosome_constant, fasta_number, five_prime_trim):
    """
    use bowtie with mm mismatches to find sequences
    """
    chromosomes = getattr(library_constants, chromosome_constant)
    if fasta_number != '1':
        for chr in chromosomes:
            subprocess.call(
                    ['bowtie', index_prefix + chr, # specify index
                         '-f', fasta_in + chr + '.fasta', # input is fasta
                         '-S', # output is sam
                         '-n', str(mm), # allow mismatches
                         '-p', '8', # parallel search threads
                         '-k', '10', # allow 10 alignments per read
                         fasta_in + chr + '.sam'])
    else:
        for chr in chromosomes:
            subprocess.call(
                    ['bowtie', index_prefix + chr, # specify index
                         '-f', fasta_in + '.fasta', # input is fasta
                         '-S', # output is sam
                         '-n', str(mm), # allow mismatches
                         '-p', '8', # parallel search threads
                         '-k', '10', # allow 10 alignments per read
                         fasta_in + chr + '.sam'])


@cli.command()
@click.argument('fasta_prefix')
@click.argument('embl_gtf')
@click.argument('nt')
@click.argument('outfile')
@click.argument('chromosome_constant')
def confirm_features(fasta_prefix, embl_gtf, nt, outfile, chromosome_constant):
    """
    find unique alignments from bwa output.
    outfile = new, multi or not_found
    nt = range of nts (a-b)
    """
    chromosomes = getattr(library_constants, chromosome_constant)
    if embl_gtf.endswith('.gtf'):
        ann_dict = parse_gtf(embl_gtf)
    elif embl_gtf.endswith('.yaml'):
        with open(embl_gtf, 'r') as file:
            contents = file.read()
        ann_dict = yaml.load(contents)
    strand = {'+':'1',
              '-':'-1',
              '0':'1',
              '16':'-1'}
    aln_dict = {}
    for chr in chromosomes:
        with open(fasta_prefix + chr + '.sam', 'r') as file:
            contents = file.readlines()
        sam_lines = [line.strip().split('\t') for line in contents \
                     if line.startswith("@") == False]
        for aln in sam_lines:
            if aln[1] == '0':
                details = {'chr':chr,
                           'start':int(aln[3]) - 1,
                           'end':int(aln[3]) + len(aln[9]),
                           'strand':ann_dict[aln[0]][4]
                           }
                if aln_dict.get(aln[0]):
                    aln_dict[aln[0]].append(details)
                else:
                    aln_dict[aln[0]] = [details]
    multi = {}
    new_anns = {}
    for feature_id in aln_dict:
        if len(aln_dict[feature_id]) > 1:
            multi[feature_id] = ann_dict[feature_id]
        elif len(aln_dict[feature_id]) == 1:
            new_anns[feature_id] = aln_dict[feature_id][0]
    if outfile == 'new':
        print yaml.dump(new_anns)
    elif outfile == 'multi':
        print yaml.dump(multi)
    elif outfile == 'not_found':
        not_found = {}
        for feature_id in ann_dict:
            if aln_dict.get(feature_id):
                continue
            else:
                if ann_dict[feature_id][0] == 'exon':
                    nt_range = nt.split('-')
                    length = ann_dict[feature_id][3] - ann_dict[feature_id][2]
                    if int(nt_range[0]) <= length:
                        if int(nt_range[1]) >= length:
                            not_found[feature_id] = ann_dict[feature_id]
        print yaml.dump(not_found)

@cli.command()
@click.argument('mapped')
@click.argument('exon_number_file')
def check_gene_representation(mapped, exon_number_file):
    """
    identify which genes are lacking representation by mapped exons
    """
    with open('/home/neil/deskgen_projects/rat_genome/Rnor_gene_list.txt', 'r') as file:
        contents = file.read()
    genes = contents.strip().split('\n')
    with open(mapped, 'r') as file:
        contents = file.read()
    all_mapped = yaml.load(contents)
    gene_rep = {}
    for gene in genes:
        n = 0
        for exon in all_mapped:
            if exon.startswith(gene):
                n += 1
        gene_rep[gene] = [n]
    with open(exon_number_file, 'r') as file:
        contents = file.readlines()
    exon_numbers = {i.strip().split('\t')[0]:[i.strip().split('\t')[1]] for i in contents}
    for gene in exon_numbers:
        exon_numbers[gene].append(str(gene_rep[gene][0]))
        exon_numbers[gene].append(str(int(exon_numbers[gene][0])-int(exon_numbers[gene][1])))
        exon_numbers[gene].append(str(float(exon_numbers[gene][1])/float(exon_numbers[gene][0])))
    lines=[]
    for i in exon_numbers:
        lines.append(i + '\t' + '\t'.join(exon_numbers[i]))
    with open('exon_representation.txt','w') as file:
        file.write('\n'.join(lines))

@cli.command()
def count_exons():
    with open('/home/neil/biodata/embl_annotations/Rattus_norvegicus.Rnor_5.0.79.gtf','r') as file:
        contents = file.readlines()
    anns = [line.strip().split('\t') for line in contents if line.startswith('#') == False]
    with open('/home/neil/deskgen_projects/rat_genome/Rnor_gene_list.txt', 'r') as file:
        contents = file.read()
    genes = contents.strip().split('\n')
    exon_no = {}
    for gene in genes:
        n = 0
        for ann in anns:
            if ann[2] == 'exon':
                if ann[8].find(gene) > -1:
                    n += 1
        exon_no[gene] = n
    lines = []
    for i in exon_no:
        lines.append(i + '\t' + str(exon_no[i]))
    with open('exon_numbers.txt', 'w') as file:
        file.write('\n'.join(lines))

@cli.command()
@click.argument('exon_rep_file')
@click.argument('mapped')
@click.argument('not_mapped')
@click.argument('list_out')
def get_exon_lists(exon_rep_file, mapped, not_mapped, list_out):
    """
    list_out = complete, unrepresented, incomplete_yes or incomplete_no
    """
    with open(exon_rep_file, 'r') as file:
        contents = file.readlines()
    exon_details = [line.strip().split('\t') for line in contents]
    with open(mapped, 'r') as file:
        contents = file.read()
    mapped_dict = yaml.load(contents)
    with open(not_mapped, 'r') as file:
        contents = file.read()
    not_mapped_dict = yaml.load(contents)
    complete = {}
    incomplete_yes = {}
    incomplete_no = {}
    unrepresented = {}
    if list_out == 'complete':
        for gene in exon_details:
            if gene[1] == gene[2]:
                for exon in mapped_dict:
                    if exon.startswith(gene[0]):
                        complete[exon] = mapped_dict[exon]
        print yaml.dump(complete)
    if list_out == 'unrepresented':
        for gene in exon_details:
            if gene[2] == '0':
                for exon in not_mapped_dict:
                    if exon.startswith(gene[0]):
                        unrepresented[exon] = not_mapped_dict[exon]
        print yaml.dump(unrepresented)
    if list_out == 'incomplete_yes':
        for gene in exon_details:
            if gene[2] != '0' and gene[1] != gene[2]:
                for exon in mapped_dict:
                    if exon.startswith(gene[0]):
                        incomplete_yes[exon] = mapped_dict[exon]
        print yaml.dump(incomplete_yes)
    if list_out == 'incomplete_no':
        for gene in exon_details:
            if gene[2] != '0' and gene[1] != gene[2]:
                for exon in not_mapped_dict:
                    if exon.startswith(gene[0]):
                        incomplete_no[exon] = not_mapped_dict[exon]
        print yaml.dump(incomplete_no)

@cli.command()
@click.argument('complete_list')
@click.argument('gene_list')
@click.argument('embl_gtf')
def check_complete_genes(complete_list, gene_list, embl_gtf):
    """
    confirm that the order of exons is correct on the same
    chromosome for complete genes
    """
    with open(gene_list, 'r') as file:
        contents = file.read()
    genes = contents.strip().split('\n')
    with open(complete_list, 'r') as file:
        contents = file.read()
    complete = yaml.load(contents)
    ann_dict = parse_gtf(embl_gtf)
    incorrect = {}
    for gene in genes:
        gene_exons = {exon_id:complete[exon_id] \
                      for exon_id in complete \
                      if exon_id.startswith(gene)}
        gtf_exons = {exon_id:ann_dict[exon_id] \
                     for exon_id in ann_dict \
                     if exon_id.startswith(gene)}
        # organise exons into transcripts
        transcripts = {}
        for exon in gene_exons:
            transcript_id = gtf_exons[exon][5]['transcript_id']
            if transcripts.get(transcript_id):
                pass
            else:
                transcripts[transcript_id] = {}
            transcripts[transcript_id][gtf_exons[exon][5]['exon_number']] = [
                            gene_exons[exon]['chr'],
                            gene_exons[exon]['start'],
                            gene_exons[exon]['end'],
                            gene_exons[exon]['strand'],
                            exon]
        # check if exons are in correct order for each transcript
        for transcript in transcripts:
            exon_nos = [int(i) for i in transcripts[transcript]]
            exon_nos.sort()
            t_strand = transcripts[transcript_id][str(exon_nos[0])][3]
            t_chr = transcripts[transcript_id][str(exon_nos[0])][0]
            if t_strand == '1':
                for exon_n in exon_nos[1:]:
                    if transcripts[transcript][str(exon_n)][0] == t_chr:
                        if transcripts[transcript][str(exon_n)][3] == t_strand:
                            if transcripts[transcript][str(exon_n -1)][2] < transcripts[transcript][str(exon_n)][1]:
                                continue
                            else:
                                incorrect[transcripts[transcript][str(exon_n)][4]] = gene_exons[transcripts[transcript][str(exon_n)][4]]
            elif t_strand == '-1':
                for exon_n in exon_nos[:len(exon_nos)-1]:
                    if transcripts[transcript][str(exon_n)][0] == t_chr:
                        if transcripts[transcript][str(exon_n)][3] == t_strand:
                            if transcripts[transcript][str(exon_n)][2] > transcripts[transcript][str(exon_n+1)][1]:
                                continue
                            else:
                                incorrect[transcripts[transcript][str(exon_n)][4]] = gene_exons[transcripts[transcript][str(exon_n)][4]]
    print yaml.dump(incorrect)

def get_complete_gene_list(complete_list, incorrect_list):
    with open(complete_list, 'r') as file:
        contents = file.read()
    complete = yaml.load(contents)
    gene_complete = set([exon_id.split('|')[0] for exon_id in complete])
    with open(incorrect_list, 'r') as file:
        contents = file.read()
    incorrect = yaml.load(contents)
    gene_incorrect = set([exon_id.split('|')[0] for exon_id in incorrect])
    complete_genes = []
    for gene in gene_complete:
        if gene in gene_incorrect:
            continue
        else:
            complete_genes.append(gene)
    with open('complete_gene_list.txt', 'w') as file:
        file.write('\n'.join(complete_genes))


@cli.command()
@click.argument('complete_gene_list')
@click.argument('complete_exons')
@click.argument('embl_gtf')
def get_cdss_from_exons(complete_gene_list, complete_exons, embl_gtf):
    """
    map cdss by comparing to exon positions,
    if both boundaries are exactly the same, they are middle exons
    otherwise they are start or end exon and can be adjusted
    """
    with open(complete_exons, 'r') as file:
        contents = file.read()
    complete = yaml.load(contents)
    ann_dict = parse_gtf(embl_gtf)
    with open(complete_gene_list, 'r') as file:
        contents = file.read()
    genes = contents.strip().split('\n')
    cdss = {}
    not_same = []
    for gene in genes:
        gene_exons = {exon_id:complete[exon_id] \
                      for exon_id in complete \
                      if exon_id.startswith(gene)}
        gtf_cdss = {exon_id:ann_dict[exon_id] \
                     for exon_id in ann_dict \
                     if exon_id.startswith(gene + '|ENSRNOP')}
        gtf_exons = {exon_id:ann_dict[exon_id] \
                     for exon_id in ann_dict \
                     if exon_id.startswith(gene + '|ENSRNOE')}
        for exon in gene_exons:
            exon_start = gtf_exons[exon][2]
            exon_end = gtf_exons[exon][3]
            transcript = gtf_exons[exon][5]['transcript_id']
            for cds in gtf_cdss:
                if transcript == gtf_cdss[cds][5]['transcript_id']:
                    cds_start = gtf_cdss[cds][2]
                    cds_end = gtf_cdss[cds][3]
                    # middle exon
                    if cds_start == exon_start and cds_end == exon_end:
                        cdss[cds] = [
                               gene_exons[exon]['chr'],
                               gene_exons[exon]['start'],
                               gene_exons[exon]['end'],
                               gene_exons[exon]['strand']]
                    # 5' exon
                    elif cds_start > exon_start and cds_end == exon_end:
                        cdss[cds] = [
                               gene_exons[exon]['chr'],
                               gene_exons[exon]['end'] - (exon_end-exon_start+1),
                               gene_exons[exon]['end'],
                               gene_exons[exon]['strand']]
                    # 3' exon
                    elif cds_start == exon_start and cds_end < exon_end:
                        cdss[cds] = [
                               gene_exons[exon]['chr'],
                               gene_exons[exon]['start'],
                               gene_exons[exon]['start'] + (exon_end-exon_start+1),
                               gene_exons[exon]['strand']]
    #with open('no_matching_CDS.txt', 'w') as file:
    #    file.write('\n'.join(not_same))
    print yaml.dump(cdss)

@cli.command()
@click.argument('incomplete_list')
@click.argument('embl_gtf')
@click.argument('biodata_prefix_fischer')
@click.argument('biodata_prefix_rnor')
@click.argument('chromosome_constant')
def align_transcripts(incomplete_list, embl_gtf, biodata_prefix_fischer, biodata_prefix_rnor, chromosome_constant):
    """
    get hypothetical transcripts to search for other exons
    """
    with open(incomplete_list, 'r') as file:
        contents = file.read()
    incomplete = yaml.load(contents)
    ann_dict = parse_gtf(embl_gtf)
    genes = set([id.split('|')[0] for id in incomplete])
    new_exon_ann = {}
    transcripts = {}
    for gene in genes:
        gene_exons = {exon_id:incomplete[exon_id] \
                      for exon_id in incomplete \
                      if exon_id.startswith(gene)}
        gtf_exons = {exon_id:ann_dict[exon_id] \
                     for exon_id in ann_dict \
                     if exon_id.startswith(gene + '|ENSRNOE')}
        gtf_transcripts = {exon_id:ann_dict[exon_id] \
                           for exon_id in ann_dict \
                           if exon_id.startswith(gene + '|ENSRNOT')}
        for exon in gene_exons:
            tranny = gene + '|' + gtf_exons[exon][5]['transcript_id']
            if transcripts.get(tranny) == None:
                start_diff = gtf_exons[exon][2] - gtf_transcripts[tranny][2] + 1000
                end_diff = gtf_transcripts[tranny][2] - gtf_exons[exon][3] + 1000
                transcripts[tranny] = {'chr':gtf_transcripts[tranny][1],
                                       'start_Fischer':gene_exons[exon]['start'] - start_diff,
                                       'end_Fischer':gene_exons[exon]['end'] + end_diff,
                                       'start_Rnor':gtf_transcripts[tranny][2],
                                       'end_Rnor':gtf_transcripts[tranny][3],
                                       'strand':gtf_transcripts[tranny][4],
                                       'exon':exon}
    chromosomes = getattr(library_constants, chromosome_constant)
    for chr in chromosomes:
        with open(biodata_prefix_fischer + chr + '.fasta') as file:
            sequence = file.read()
        chr_fasta_Fischer = sequence.strip().split('\n')[1]
        with open(biodata_prefix_rnor + chr + '.fasta') as file:
            sequence = file.read()
        chr_fasta_Rnor = sequence.strip().split('\n')[1]
        chr_transcripts = {t:transcripts[t] for t in transcripts \
                      if transcripts[t]['chr'] == chr}
        for trans in chr_transcripts:
            with open('transcript_sequences_Fischer.fasta', 'w') as file:
                file.write('>' + '|'.join([trans, 
                                           chr, 
                                           str(transcripts[trans]['start_Fischer']),
                                           str(transcripts[trans]['end_Fischer']),
                                           transcripts[trans]['strand'],
                                           transcripts[trans]['exon'],
                                           'Fischer']) + '\n')
                file.write(chr_fasta_Fischer[transcripts[trans]['start_Fischer']:
                                             transcripts[trans]['end_Fischer']] + '\n')
            with open('transcript_sequences_Rnor.fasta', 'w') as file:
                file.write('>' + '|'.join([trans, 
                                           chr, 
                                           str(transcripts[trans]['start_Rnor']),
                                           str(transcripts[trans]['end_Rnor']),
                                           transcripts[trans]['strand'],
                                           transcripts[trans]['exon'],
                                           'Rnor']) + '\n')
                file.write(chr_fasta_Rnor[transcripts[trans]['start_Rnor']-1:
                                             transcripts[trans]['end_Rnor']] + '\n')
            subprocess.call(['needle', 
                             '-asequence', 'transcript_sequences_Fischer.fasta',
                             '-bsequence', 'transcript_sequences_Rnor.fasta',
                             '-gapopen', '1.0', 
                             '-gapextend', '0.5',
                             '-outfile', 'Fischer-Rnor.align'])
            alignment = AlignIO.read('Fischer-Rnor.align', 'emboss')
            fischer_seq = str(alignment.get_all_seqs()[0].seq)
            rnor_seq = str(alignment.get_all_seqs()[1].seq)
            fischer_pos = transcripts[trans]['start_Fischer']
            rnor_pos = transcripts[trans]['start_Rnor']-1
            gap_dict = {rnor_pos:fischer_pos}
            for nt in range(1,len(fischer_seq),1):
                if fischer_seq[nt] != '-':
                    fischer_pos += 1
                if rnor_seq[nt] != '-':
                    rnor_pos += 1
                if gap_dict.get(rnor_pos) == None:
                    gap_dict[rnor_pos] = [fischer_pos]
            for exon in gtf_exons:
                if trans.split('|')[0] == gtf_exons[exon][5]['transcript_id']:
                    new_exon_ann[exon] = {'chr':gtf_exons[exon][1],
                                          'start':gap_dict[gtf_exons[exon][2]],
                                          'end':gap_dict[gtf_exons[exon][3]],
                                          'strand':gtf_exons[exon][4]
                                          }
    print yaml.dump(new_exon_ann)


if __name__ == '__main__':
    cli()