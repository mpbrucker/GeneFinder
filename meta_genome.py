from load import load_nitrogenase_seq, load_metagenome
import doctest

nitrogenase = load_nitrogenase_seq()


def longest_common_substring(string1, string2):
    """
    Computes the longest common substring between two strings.

    >>> longest_common_substring("qwerrasldjfasdf","palljhasldjfpljhplk")
    ['asldjf']

    >>> longest_common_substring("oifusadkdbd","asdfoifusadkdbdwrert")
    ['oifusadkdbd']

    >>> longest_common_substring("test","test1")
    ['test']

    >>> longest_common_substring("abc","def")
    []

    >>> longest_common_substring("xabcxdef","jabcjdef")
    ['abc', 'def']
    """
    substring_array = [[0 for x in range(len(string1))] for y in range(len(string2))]
    max_sub = 0  # length of longest substring
    longest = []  # The array of longest substrings
    for y in range(0, len(string2)):
        for x in range(0, len(string1)):
            if string1[x] == string2[y]:  # if the current letter in the two strings match
                if x == 0 or y == 0:
                    substring_array[y][x] = 1
                else:
                    substring_array[y][x] = substring_array[y-1][x-1] + 1
                if substring_array[y][x] > max_sub:  # There is a new longest substring
                    max_sub = substring_array[y][x]
                    longest = [string1[x-max_sub+1:x+1]]
                elif substring_array[y][x] == max_sub:  # The substring is of equal length
                    longest.append(string1[x-max_sub+1:x+1])
            else:
                substring_array[y][x] == 0
    return longest


def find_possible_nitrogenase():
    """
    Finds the best match for the nitrogenase genome out of the metagenome by finding
    the genome with the longest common substring between the two.
    """
    max_len = 0
    metagenome = load_metagenome()
    # Nice to use the line below for demonstration purposes, otherwise it takes a LONG time
    # metagenome = metagenome[:25]
    for genome in metagenome:  # Iterate through all genomes in the metagenome
        lc_substring = longest_common_substring(nitrogenase, genome[1])  # Get longest substring
        if len(lc_substring) > max_len:  # there is a new longest common substring
            max_len = len(lc_substring)
            best_match = genome  # There is a new best match
    return best_match


if __name__ == '__main__':
    doctest.testmod()
    best = find_possible_nitrogenase()
    print(best)
