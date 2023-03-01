from constants import *
from numpy import array, delete, where


def needleman_wunsch(seq1, seq2, matrix, gap=-1) -> int:
    """
    Algorithme de Needleman-Wunsch pour aligner deux séquences

    Arguments:
    seq1 -- première séquence (str)
    seq2 -- deuxième séquence (str)
    gap -- score de trou (int, default=-1)

    Retourne:
    score -- score d'alignement (int)
    """

    # Initialisation de la matrice de scores
    n = len(seq1)
    m = len(seq2)
    score_matrix = [[0] * (m+1) for i in range(n+1)]

    # Initialisation de la première ligne et colonne de la matrice
    for i in range(1, n+1):
        score_matrix[i][0] = gap * i
    for j in range(1, m+1):
        score_matrix[0][j] = gap * j

    # Remplissage de la matrice de scores
    for i in range(1, n+1):
        for j in range(1, m+1):
            indice_lig = DAA[seq1[i-1]]
            indice_col = DAA[seq2[j-1]]
            match_mismatch_score = matrix[indice_lig, indice_col]
            score_matrix[i][j] = max(
                score_matrix[i-1][j-1] +
                match_mismatch_score,  # Correspondance
                score_matrix[i-1][j] + gap,  # Trou dans seq1
                score_matrix[i][j-1] + gap  # Trou dans seq2
            )

    return score_matrix[n][m]


def open_fasta(filename) -> dict:
    """open fasta and return all the sequence in a dictionnary

    Args:
        filename (str): path to file

    Returns:
        dict: all the seq of the file
    """
    list_of_seq = {}
    acid_amine = ['C', 'S', 'T', 'A', 'G', 'P', 'D', 'E', 'Q',
                  'N', 'H', 'R', 'K', 'M', 'I', 'L', 'V', 'W', 'Y', 'F']
    with open(filename, 'r') as file:
        for line in file:
            if not line.startswith(">") and line.strip().split() in acid_amine:
                return False
            elif line.startswith(">"):
                temp_header = line.strip()
                list_of_seq[line.strip()] = ""
            else:
                list_of_seq[temp_header] += line.strip()
        return list_of_seq


def get_all_max_score(filename):
    """get the max score of all sequence

    Args:
        filename (str): path of the file

    Returns:
        list: array of array with all score
    """
    counter_col = 0
    counter_lig = 0
    list_name_clades = []
    fasta_dict = open_fasta(filename)
    index_dict = len(fasta_dict)
    matrix_distance = array([[0] * (index_dict)
                            for _ in range(index_dict+1)])
    for i in range(len(matrix_distance)+1):
        list_name_clades.append("seq"+str(i+1))

    for header in fasta_dict.keys():
        counter_col = 0
        counter_lig += 1
        seq_one = fasta_dict[header]
        for key, seq_two in fasta_dict.items():
            if counter_col >= counter_lig:
                break
            if header != key:
                matrix_distance[counter_lig, counter_col] = needleman_wunsch(
                    seq_one, seq_two, BLOSUM62)
            counter_col += 1
    matrix_distance = delete(matrix_distance, 0, 0)
    matrix_distance = delete(matrix_distance, 0, 0)

    return matrix_distance


def find_closer_upgma(matrix_distance):
    """find the two clades to merge first by finding
    the minimal distance in the given matrix.

    Args:
        matrix_distance (list): Tableau des distances

    Returns:
        int: the minimum of the all the score
    """
    mini = 100000
    for x in range(0, len(matrix_distance)):
        for y in range(0, len(matrix_distance)):
            if x != y:
                temp_min = matrix_distance[x, y]
                if temp_min < mini and temp_min != 0 :
                    mini, line_min, col_min = temp_min, x+1, y+1
    return mini, line_min, col_min


if __name__ == "__main__":
    print("\n", get_all_max_score("opsines_juste4.fasta.txt"))
    print("\n", find_closer_upgma(get_all_max_score("opsines_juste4.fasta.txt")))
