import pytest

@pytest.fixture
def fasta_sequence(request):
    """
	Make a mock fasta file of chr input
	"""
    return 'ATGCTAGCTAGCTAGCTAGCTAGCGATGCTAGCTGCATGCTAGCTA'\
    'GCATGCAGCTGCATCGATCGATCGATGCATGCTAGCTAGATGTCGATCGATGCA'

@pytest.fixture
def snp_table(request):
    """
    A mock input file of snps
    """
    f = open('snp_test.vcf', 'w')
    f.write(
    '1\t10\t1\tA\tT\n1\t25\t2\tG\tC\n1\t79\t3\tT\tC\n'
    )
    f.close

@pytest.fixture
def indel_table(request):
    """
    A mock input file of indels
    """
    f = open('indel_test.vcf', 'w')
    f.write(
    '1\t10\t1\tA\tTCTG\n1\t25\t2\tGATGCTAG\tC\n1\t79\t3\tT\tCCCTGTAAA\n'
    )

