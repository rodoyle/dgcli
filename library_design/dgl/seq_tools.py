import click

@click.group()
def cli():
    pass

CODON_TABLE = {
               'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
               'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
               'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
               'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
               'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
               'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
               'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
               'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
               'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
               'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
               'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
               'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
               'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
               'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
               'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
               'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
               }

COMPLEMENT = {
              'A':'T',
              'T':'A',
              'C':'G',
              'G':'C',
              'N':'N'
              }

@cli.command()
@click.argument('na_seq')
@click.argument('all', default=False)
def translate(na_seq, all):
    """
    clickable function
    Return amino acid sequence string of nucleic acid sequence (na_seq) input using
    phase 0. na_seq can be DNA (ACGT) or RNA (ACGU)
    if all == True, will return a list of 6 sequences - translation of all 
    forward and reverse phases [0F, 1F, 2F, 0R, 1R, 2R]
    if no translation for codon, script will return '-'
    """
    if na_seq.find('U') > -1:
        na_seq = na_seq.replace('U','T')
    if all == False:
        aa = ''
        for i in range(0, len(na_seq), 3):
            if CODON_TABLE.get(na_seq[i:i+3]):
                aa += CODON_TABLE[na_seq[i:i+3]]
            else:
                aa += '-'
        print aa
    else:
        seq_out = []
        for j in range(3):
            aa = ''
            for i in range(j, len(na_seq), 3):
                if CODON_TABLE.get(na_seq[i:i+3]):
                    aa += CODON_TABLE[na_seq[i:i+3]]
                else:
                    aa += '-'
            seq_out.append(aa)
        rev_seq = ''.join([COMPLEMENT[i] for i in reversed(na_seq)])
        for j in range(3):
            aa = ''
            for i in range(j, len(rev_seq), 3):
                if CODON_TABLE.get(rev_seq[i:i+3]):
                    aa += CODON_TABLE[rev_seq[i:i+3]]
                else:
                    aa += '-'
            seq_out.append(aa)
        print seq_out

def translate_function(na_seq, all=False):
    """
    non-clickable, callable function
    Return amino acid sequence string of nucleic acid sequence (na_seq) input using
    phase 0. na_seq can be DNA (ACGT) or RNA (ACGU)
    if all == True, will return a list of 6 sequences - translation of all 
    forward and reverse phases [0F, 1F, 2F, 0R, 1R, 2R]
    if no translation for codon, script will return '-'
    """
    if na_seq.find('U') > -1:
        na_seq = na_seq.replace('U','T')
    if all == False:
        aa = ''
        for i in range(0, len(na_seq), 3):
            if CODON_TABLE.get(na_seq[i:i+3]):
                aa += CODON_TABLE[na_seq[i:i+3]]
            else:
                aa += '-'
        return aa
    else:
        seq_out = []
        for j in range(3):
            aa = ''
            for i in range(j, len(na_seq), 3):
                if CODON_TABLE.get(na_seq[i:i+3]):
                    aa += CODON_TABLE[na_seq[i:i+3]]
                else:
                    aa += '-'
            seq_out.append(aa)
        rev_seq = ''.join([COMPLEMENT[i] for i in reversed(na_seq)])
        for j in range(3):
            aa = ''
            for i in range(j, len(rev_seq), 3):
                if CODON_TABLE.get(rev_seq[i:i+3]):
                    aa += CODON_TABLE[rev_seq[i:i+3]]
                else:
                    aa += '-'
            seq_out.append(aa)
        return seq_out

if __name__ == '__main__':
    cli()