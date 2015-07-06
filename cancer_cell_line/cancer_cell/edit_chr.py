import subprocess

from Bio import SeqIO

def fasta_single_seq_parser(fasta_file):
    """
    return a single fasta sequence as a string for a single-sequence fasta file
    e.g. a chromosome sequence
    :param fasta_file: fasta file name
    :return:
    """
    if fasta_file.endswith('.gz'):
        subprocess.call(
            ['gzip', '-d', '-k', fasta_file]
        )
        fasta_name = fasta_file.split('.gz')[0]
    else:
        fasta_name = fasta_file
    with open(fasta_name) as file:
        fasta_lines = file.readlines()
    fasta_seq = ''
    for line in fasta_lines[1:]:
        fasta_seq += line.strip()
    return fasta_seq


def add_snps(fasta_file, chr, genome, snp_table):
    """
    use caveman.vcf-type file to apply SNPs to
    chromosome sequence
    """
    fasta = fasta_single_seq_parser(fasta_file)
    dodgy = []
    with open(snp_table) as file:
        contents = file.readlines()
    snps = [line.strip().split('\t') for line in contents if line.startswith(chr + '\t') == True]
    for snp in snps:
        if int(snp[1]) <= len(fasta):
            if fasta[int(snp[1])-1] == snp[3]:
                fasta = fasta[0:int(snp[1])-1] + snp[4] + fasta[int(snp[1]):]
            else:
                dodgy.append(
                    'incorrect reference sequence: ' + snp[0] + ' ' + snp[1] + ' ' + snp[2] + ' ' + fasta[int(snp[1])-1])
    fasta_file = open('chr' + chr + '_' + genome + '_snps.flat', 'w')
    #fasta_file.write('>SNPed_chr' + chr + '_' + genome + '\n' + fasta + '\n')
    fasta_file.write(fasta)
    fasta_file.close()
    if len(dodgy) > 0:
        error_file = open('chr' + chr + '_' + genome + '_incorrect_snps.txt', 'w')
        error_file.write('\n'.join(dodgy))
        error_file.close()

def map_indels(chr_fasta, chr, genome, indel_table, outfile_prefix):
    """
    use pindel.vcf-type file to integrate INDELs into
    chromosome sequence and make a  cumulative map file
    to adjust feature coordinates.
    """
    fasta = str(SeqIO.read(chr_fasta, 'fasta').seq)
    with open(indel_table) as file:
        contents = file.readlines()
    indels = [line.strip().split('\t') for line in contents if line.startswith(chr + '\t') == True]
    map = []
    dodgy = []
    # start from end of chromosome to retain coordinates
    for indel in reversed(indels):
        if int(indel[1]) <= len(fasta):
            if fasta[int(indel[1])-1] == indel[3][0]:
                fasta = fasta[0:int(indel[1])-1] + indel[4] + fasta[(int(indel[1])-1) + len(indel[3]):]
                map.append([indel[0], indel[1], str(len(indel[4]) - len(indel[3]))])
            else:
                dodgy.append(
                'incorrect reference sequence: ' + indel[0] + ' ' + indel[1] + ' ' + indel[2] + ' ' + fasta[int(indel[1])-1])
    fasta_file = open(outfile_prefix + '.fasta', 'w')
    fasta_file.write('>INDELed_chr' + chr + '\n' + fasta + '\n')
    fasta_file.close()
    if len(dodgy) > 0:
        error_file = open('chr' + chr + '_' + genome + '_incorrect_INDELs.txt', 'w')
        error_file.write('\n'.join(dodgy))
        error_file.close()
    map_cum = []
    for change in reversed(map):
        if len(map_cum) == 0:
            map_cum.append(change)
        else:
            map_cum.append([change[0], change[1], str(int(map_cum[len(map_cum)-1][2]) + int(change[2]))])
    map_line = ['\t'.join(entry) for entry in map_cum]
    map_file = open(outfile_prefix + '.map', 'w')
    map_file.write('\n'.join(map_line))
    map_file.close()

def for_all_chr():
    chromosomes = os.listdir('.')
    for n,chr in enumerate(chromosomes[1:]):
        #subprocess.call(['gzip', '-d', chr])
        #chr_unzip = chr.split('.gz')
        chr_n = chr.split('chr')
        chr_no = chr_n[1].split('.')
        add_snps(chr, chr_no[0], 'GRCh38', '../A375_GRCh38.p2_remap_caveman.vcf')



# chr_fasta = 'Homo_sapiens.GRCh37.75.dna_rm.chromosome.20.fa'
# snp_table = 'A375_caveman.vcf'
# chr_fasta = 'chr20_SNPed.fasta'
# indel_table = 'A375_pindel.vcf'

# edit_chr.add_snps('Homo_sapiens.GRCh37.75.dna_rm.chromosome.20.fa', '20', 'GRCh37', 'A375_caveman.vcf', 'chr20_GRCh37_SNPed.fasta')
# edit_chr.map_indels('chr20_GRCh37_SNPed.fasta', '20', 'GRCh37', 'A375_pindel.vcf', 'chr20_GRCh37_INDELs')

# edit_chr.add_snps('Homo_sapiens.GRCh38.dna_rm.chromosome.20.fa', '20', 'GRCh38', 'A375_caveman.vcf', 'chr20_GRCh38_SNPed.fasta')
# edit_chr.map_indels('chr20_GRCh38_SNPed.fasta', '20', 'GRCh38', 'A375_pindel.vcf', 'chr20_GRCh38_INDELs')


