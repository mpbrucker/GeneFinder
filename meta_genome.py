from load import load_nitrogenase_seq, load_metagenome
import doctest
from multiprocessing import Pool, Lock, Manager
from functools import partial


def longest_common_substring(string2, lock, old_longest, in_string):
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
    # print('Nitrogenase: ',string2)
    string1 = in_string[1]
    # print('Metagenome: ',string1)
    # Builds blank array of dimensions str1 len x str 2 len
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
    lock.acquire()
    # print(longest)
    if (len(longest)) > len(old_longest[1]):
        old_longest[0] = in_string[0]
        old_longest[1] = longest
    # print(old_longest)
    lock.release()


def find_possible_nitrogenase():
    """
    Finds the best match for the nitrogenase genome out of the metagenome by finding
    the genome with the longest common substring between the two.
    """
    nitrogenase = load_nitrogenase_seq()
    metagenome = load_metagenome()
    print("Length: ", len(metagenome))
    process_pool = Pool(15)
    m = Manager()
    l = m.Lock()
    metagenome = metagenome[:25]
    cur_longest = ['', '']
    lc_sub = partial(longest_common_substring, nitrogenase, l, cur_longest)
    process_pool.map(lc_sub, metagenome)
    process_pool.close()
    process_pool.join()
    print(cur_longest)

    # for genome in metagenome:  # Iterate through all genomes in the metagenome
    #
    #     lc_substring = longest_common_substring(nitrogenase, genome[1])
    #     if len(lc_substring) > max_len:  # there is a new longest common substring
    #         max_len = len(lc_substring)
    #         best_match = genome


if __name__ == '__main__':
    find_possible_nitrogenase()
    # doctest.testmod()
