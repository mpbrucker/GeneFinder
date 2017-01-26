# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Matt Brucker

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq
import re
from itertools import repeat


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    nucleotides = ["G", "A", "T", "C"]
    # Returns the complementary corresponding nucleotide
    return nucleotides[3-nucleotides.index(nucleotide.upper())]


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    return ''.join([get_complement(val) for val in dna[::-1]])


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGCATGAATGTAGATAGTAGTGCCC")
    'ATGCATGAATGTAGA'
    """
    return dna[0:min(list(filter(lambda y: y % 3 == 0, [item for sublist in list(map(lambda x: [m.start() for m in re.finditer(x, dna)], ["TAG", "TAA", "TGA"])) for item in sublist]))+list([len(dna)]))]


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGCATATGTGTAGATAGATGTGCCC")
    ['ATGCATATGTGTAGA', 'ATGTGCCC']
    """
    end_ORF_index = 0
    current_frame = dna
    all_ORFs = []
    while len(current_frame) > 0:
        next_ORF_index = min(list(filter(lambda y: y % 3 == 0, [m.start() for m in re.finditer('ATG', current_frame)])))
        current_frame = current_frame[next_ORF_index:]
        cur_ORF = rest_of_ORF(current_frame)
        all_ORFs.append(cur_ORF)
        current_frame = current_frame[len(cur_ORF)+3:]

    return all_ORFs


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    all_ORFs = []
    for counter in range(0, 3):
        current_frame = dna[counter:]
        while len(current_frame) > 0:
            all_ORF_index = [y for y in [m.start() for m in re.finditer('ATG', current_frame)] if y % 3 == 0]
            if len(all_ORF_index) > 0:
                next_ORF_index = min(all_ORF_index)
            else:
                break
            current_frame = current_frame[next_ORF_index:]
            cur_ORF = rest_of_ORF(current_frame)
            all_ORFs.append(cur_ORF)
            current_frame = current_frame[len(cur_ORF)+3:]

    return all_ORFs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    reverse_complement = get_reverse_complement(dna)
    all_ORFs = find_all_ORFs(dna) + find_all_ORFs(reverse_complement)
    return all_ORFs


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    return max(find_all_ORFs_both_strands(dna))


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    cur_max = 0
    for idx in range(0, num_trials):
        new_max = len(longest_ORF(shuffle_string(dna)))
        if new_max > cur_max:
            cur_max = new_max

    return new_max
    # Top secret 1-line version
    # return len(max([longest_ORF(shuffle_string(val)) for val in list(repeat(dna, num_trials))]))


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    return ''.join([aa_table[codon] for codon in [dna[idx*3:idx*3+3] for idx in range(0, len(dna)//3)]])


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    print(threshold)
    all_ORF = [orf for orf in find_all_ORFs_both_strands(dna) if len(orf) > threshold]
    return [coding_strand_to_AA(val) for val in all_ORF]


if __name__ == "__main__":
    # import doctest
    # doctest.testmod()
    dna_test = load_seq("./data/X73525.fa")
    print(gene_finder(dna_test))
