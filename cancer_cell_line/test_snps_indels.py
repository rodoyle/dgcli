import pytest
import subprocess

from cancer_cell import edit_chr


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