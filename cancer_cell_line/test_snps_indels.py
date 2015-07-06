#!/usr/bin/python

import pytest

from cancer_cell.edit_chr import fasta_single_seq_parser

def test_add_snps():
	edit_chr.add_snps('chr_test.fasta', '1', 'unicorn', 'snp_test.txt', 'snp_test.fasta')
	test_fasta = str(SeqIO.read('snp_test.fasta', 'fasta').seq)
	assert test_fasta[9] == 'T'
	assert test_fasta[24] == 'C'
	assert test_fasta[78] == 'C'

def test_map_indels():
	edit_chr.map_indels('chr_test.fasta', '1', 'unicorn', 'indel_test.txt', 'indel_test')
	test_fasta = str(SeqIO.read('indel_test.fasta', 'fasta').seq)
	with open('indel_test.map', 'r') as file:
		contents = file.readlines()
	test_map = [line.strip().split('\t') for line in contents]
	coor_diff = 1
	for entry in test_map:
		coor_diff = int(entry[2]) + coor_diff
	assert test_fasta[9] == 'T'
	assert len(test_fasta) == coor_diff + 100