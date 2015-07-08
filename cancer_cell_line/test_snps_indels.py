import pytest
import subprocess
import click

#from cancer_cell import edit_chr


def test_add_snps(fasta_sequence, snp_table):
    snp_table
    fasta_out = edit_chr.add_snps(fasta_sequence, '1', 'snp_test.vcf')
    assert fasta_out[9] == 'T'
    assert fasta_out[24] == 'C'
    assert fasta_out[78] == 'C'
    assert len(fasta_out) == 100
    subprocess.call(['rm', 'snp_test.vcf'])

def test_map_indels(fasta_sequence, indel_table):
    indel_table
    fasta_out = edit_chr.add_indels_with_mapping(fasta_sequence, '1', 'indel_test.vcf')
    with open('indel_map.txt', 'r') as file:
        contents = file.readlines()
    test_map = [line.strip().split('\t') for line in contents]
    coor_diff = 1
    for entry in test_map:
        coor_diff = int(entry[2]) + coor_diff
    assert fasta_out[9] == 'T'
    assert len(fasta_out) == coor_diff + 100
    subprocess.call(['rm', 'indel_map.txt'])
    subprocess.call(['rm', 'indel_test.vcf'])

@click.group()
def cli():
    pass

@cli.command()
@click.argument('vcf_file')
def check_for_snp_indel_collisions(vcf_file):
    with open(vcf_file) as file:
        contents = file.readlines()
    vcf = [i.strip().split('\t') for i in contents if i.startswith('#') == False]
    colo = ['chromosome_name,start,end,strand,label']
    for v in vcf:
        colo.append(','.join(
            [v[0], # chromosome
             v[1], # start
             str(int(v[1])+len(v[3])), # end
             '1', # strand
             '|'.join(v[0:3]), # header
             ]
        ))
    with open('colo_vcf.csv', 'w') as file:
        file.write('\n'.join(colo) + '\n')
    outfile = open('colo_out.csv', 'w')
    subprocess.check_call(
        ['colocator',
         '--search', 'colo_vcf.csv',
         '--target', 'colo_vcf.csv',
        ],
        stdout=outfile)
    outfile.close()
    with open('colo_out.csv') as file:
        contents = file.readlines()
    colocated = [i.strip().split(',') for i in contents]
    duplicated = []
    for colo in colocated:
        if colo[0] != colo[1]:
            duplicated.append(','.join(colo))
    subprocess.call(['rm', 'colo_vcf.csv'])
    subprocess.call(['rm', 'colo_out.csv'])
    if len(duplicated) == 0:
        print 'No Collisions Detected'
    else:
        print '\n'.join(duplicated)


if __name__ == '__main__':
    cli()
