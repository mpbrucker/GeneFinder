# -*- coding: utf-8 -*-
"""
Analyzes strings of nucleotides to find chains of amino acids.
Based on an imported string of nucleotides, finds the most likely
coded amino acids.

@author: Matt Brucker

"""

import random
from amino_acids import codons, aa_table   # you may find these useful
from load import load_seq
import re


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
    # The list of nucleotide.  Each complement is on the opposite end of the list as the original nucleotide.
    nucleotides = ["G", "A", "T", "C"]
    # Returns the same nucleotide except indexed backward (i.e. the complement)
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
    # Reverses dna and builds a list of the complement of each letter, then joining it into a string.
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
    # Builds a list of all indices of the stop codons in the string
    all_end_indices = list(map(lambda x: [m.start() for m in re.finditer(x, dna)], codons[10]))
    # Flattens the list
    flat_end_indices = [item for sublist in all_end_indices for item in sublist]
    # Filters out indices that aren't valid (i.e. not divisible by 3)
    valid_end_indices = list(filter(lambda y: y % 3 == 0, flat_end_indices))
    # Adds the end index (in case there are no stop codons and frame) and finds the minimum
    end_index = min(valid_end_indices+list([len(dna)]))
    # Returns the substring of dna from the beginning to the stopping point
    return dna[0:end_index]


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
    current_frame = dna
    all_ORFs = []
    while len(current_frame) > 0:
        # Finds the beginning index of the next valid ORF
        all_possible_ORF = [ind.start() for ind in re.finditer(codons[3][0], current_frame)]
        valid_indices = list(filter(lambda val: val % 3 == 0, all_possible_ORF))
        if len(valid_indices) == 0:  # There are no valid ORFs left in frame
            break
        else:
            next_ORF_index = min(valid_indices)
        current_frame = current_frame[next_ORF_index:]  # slices the frame to start at the ORF

        # Adds the ORF to a list and slices the original frame to start at the end of the ORF
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
        current_frame = dna[counter:]  # Shifts the string to start at the shifted frame
        all_ORFs.extend(find_all_ORFs_oneframe(current_frame))

    return all_ORFs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # Concats the list of all ORFs for the original and reverse complement strands
    reverse_complement = get_reverse_complement(dna)
    all_ORFs = find_all_ORFs(dna) + find_all_ORFs(reverse_complement)
    return all_ORFs


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # Returns the longest ORF.  If two ORFs are the same length, it uses the numerical values of characters
    return max(find_all_ORFs_both_strands(dna))


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """

    cur_max = 0
    for idx in range(0, num_trials):
        # Find the longest ORF in this random shuffle
        new_max = len(longest_ORF(shuffle_string(dna)))
        if new_max > cur_max:
            cur_max = new_max  # A new max has been found, set it

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
    # Finds the indices of all codons, and generates a list of the three-letter substrings starting at those indices
    codons = [dna[idx*3:idx*3+3] for idx in range(0, len(dna)//3)]
    # Returns a list of the amino acids corresponding to each codon
    return ''.join([aa_table[codon] for codon in codons])


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500)  # Get a reasonable number for minimum length
    # Gets all ORFs with length over the threshold
    all_ORF = [orf for orf in find_all_ORFs_both_strands(dna) if len(orf) > threshold]
    # Returns a list of the AAs corresponding to the ORFs in all_ORF
    return [coding_strand_to_AA(val) for val in all_ORF]


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
    dna = load_seq("./data/X73525.fa")
    print(gene_finder(dna))
