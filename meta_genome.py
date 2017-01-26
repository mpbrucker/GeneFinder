from load import load_nitrogenase_seq, load_metagenome
import doctest


def longest_common_substring(string1, string2):
    """
    Computes the longest common substring between two strings.

    >>> longest_common_substring("qwerrasldjfasdf","palljhasldjfpljhplk")
    ["asldjf"]

    >>> longest_common_substring("oifusadkdbd","asdfoifusadkdbdwrert")
    ["oifusadkdbd"]
    """
    substring_array = [[0 for x in range(len(string1))] for y in range(len(string2))]
    max_sub = 0
    longest = []
    # try:
    for i in range(0, len(string2)):
        for j in range(0, len(string1)):
            if string1[j] == string2[i]:
                if i == 0 or j == 0:
                    substring_array[i][j] == 1
                else:
                    substring_array[i][j] = substring_array[i-1][j-1] + 1

                if substring_array[i][j] > max_sub:
                    max_sub = substring_array[i][j]
                    print('Max sub:',max_sub)
                    print('x:',i)
                    print('y:',j)
                    longest = [string1[i-max_sub+1:i]]
                # elif substring_array[i][j] == max_sub:
                #     longest.append(string1[i-max_sub+1:i])
            else:
                substring_array[i][j] == 0
    for row in substring_array:
        print(row)
    return longest
    # except IndexError:
    #     print(i)
    #     print(j)


if __name__ == '__main__':
    print('Test')
    doctest.testmod()
