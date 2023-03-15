from constants import *
import numpy as np
#from upgma import Upgma

# def freq_matrix(seq_list):
#     """Cette fonction prend une liste de séquences d'acides aminés et retourne une matrice de fréquence"""
#     num_seq = len(seq_list)
#     seq_len = len(seq_list[0])
#     freq_matrix = np.zeros((seq_len, 20))
#     for i in range(seq_len):
#
   #for j in range(num_seq):
#             if seq_list[j][i] != "-":
#                 freq_matrix[i][ord(seq_list[j][i])-65] += 1
#     freq_matrix /= num_seq
#     return freq_matrix

def freq_matrix(seq_list):
    """Cette fonction prend une liste
    de séquences d'acides aminés et
    retourne une matrice de fréquence"""
    num_seq = len(seq_list)
    seq_len = len(seq_list[0])
    freq_matrix = np.zeros((seq_len, 21))
    # Liste des acides aminés
    amino_acids = "CSTAGPDEQNHRKMILVWYF-"
    for i in range(seq_len):
        for j in range(num_seq):
            indice_aa_matrice = amino_acids.find(seq_list[j][i])
            freq_matrix[i][indice_aa_matrice] += 1
    freq_matrix /= num_seq
    return freq_matrix


def needleman_wunsch_profile(seq1 : list, seq2 : list, matrix : list, gap=-4):
    """Cette fonction prend en entrée deux profils
    de séquences multiples, une matrice de score,
    et une valeur de gap (par défaut -4) et
    renvoie un nouveau profil aligné"""

    if type(seq1) is not list:
        seq1 = [seq1]
    if  type(seq2) is not list:
        seq2 = [seq2]

    m = len(seq1[0]) + 1
    n = len(seq2[0]) + 1

    freq_matrix1 = freq_matrix(seq1)
    freq_matrix2 = freq_matrix(seq2)

    # Initialisation de la matrice des scores et de direction
    score_matrix = np.zeros((m, n))
    direction_matrix = np.zeros((m, n))
    for i in range(1, m):
        score_matrix[i][0] = score_matrix[i-1][0] + gap
        direction_matrix[i][0] = 1
    for j in range(1, n):
        score_matrix[0][j] = score_matrix[0][j-1] + gap
        direction_matrix[0][j] = 2

    print(score_matrix)
    
    # Remplissage de la matrice des scores et de direction
    for i in range(1, m):
        for j in range(1, n):
            # Calcul du score pour un alignement diagonal
            diag_score = score_matrix[i-1][j-1] + np.sum(freq_matrix1[i-1] * matrix * freq_matrix2[j-1])
            # Calcul du score pour un alignement vertical
            up_score = score_matrix[i-1][j] + gap
            # Calcul du score pour un alignement horizontal
            left_score = score_matrix[i][j-1] + gap
            # Choix du score maximum et de la direction correspondante
            if diag_score == left_score and diag_score == up_score:
                score_matrix[i][j], direction_matrix[i][j] = (diag_score, 0)

            elif diag_score == left_score and diag_score != up_score:
                score_matrix[i][j], direction_matrix[i][j] = max((diag_score, 0), (up_score, 1))

            elif diag_score != left_score and diag_score == up_score:
                score_matrix[i][j], direction_matrix[i][j] = max((diag_score, 0), (left_score, 2))

            elif diag_score != left_score and diag_score != up_score:
                score_matrix[i][j], direction_matrix[i][j] = max((diag_score, 0), (up_score, 1), (left_score, 2))
    print(score_matrix)
    print(direction_matrix)
    # Construction du nouvel alignement
    aligned_seq1 = []
    aligned_seq2 = []

    print("\n se1 : ",seq1)
    for k in seq1:
        i = m - 1
        j = n - 1
        #print('\n new_aa_to_join: ',new_aa_to_join)
        print("longueur seq1: ", len(seq1))
        new_aa_to_join = ""
        while i > 0 or j > 0:

            direction = direction_matrix[i][j]
            if direction == 0:
                new_letter = k[i-1]
                new_aa_to_join += new_letter
                i -= 1
                j -= 1

            elif direction == 1:
                new_letter = k[i-1]
                new_aa_to_join += new_letter
                i -= 1
            elif direction == 2:
                new_aa_to_join += "-"
                j -= 1
        print("\n APRES BOUCLE FOR",new_aa_to_join)
        aligned_seq1.append(new_aa_to_join[::-1])

    for k in seq2:
        i = m - 1
        j = n - 1
        new_aa_to_join2 = ""
        while i > 0 or j > 0:
            direction = direction_matrix[i][j]
            if direction == 0:
                new_letter = k[j-1]
                new_aa_to_join2 += new_letter
                i -= 1
                j -= 1

            elif direction == 1:
                new_aa_to_join2 += "-"
                i -= 1

            elif direction == 2:
                new_letter = k[j-1]
                new_aa_to_join2 += new_letter
                j -= 1
        aligned_seq2.append(new_aa_to_join2[::-1])


    print("\n aligned_seq1: ",aligned_seq1)
    new_alignment = []
    for elem in aligned_seq1:
        new_alignment.append(elem)
    for elem in aligned_seq2:
        new_alignment.append(elem)

    print("\n NEW ALIGNEMENTS: ", new_alignment)

    return new_alignment


if __name__ == "__main__":
    #p1 = ["MSKMSEE", "MSKMSEE"]
    p1 = ["MSKMSEE", "TGWSLAF"]
    # p1 = ["RQPLNYILVNVSLGGFIYCIFSVFIVFITSCYGYFVFGRHVCALEAFLGCTAGL", "VTGWSLAFLAFERYIIICKPFGNFRFSSKHALMVVVATWTIGIGVSIPPFFGGG"]
    p2 = ["MRKMSEEEFYLFKNISSVGPWDGPQYHIAPVWAFYLQAAFMGTVFLIGFPLNAMVLVATL", "RYKKLRQPLNYILVNVSFGGFLLCIFSVFPVFVASCNGYFVFGRHVCALEGFLGTVAGLV",
          "TGWSLAFLAFERYIVICKPFGNFRFSSKHALTVVLATWTIGIGVSIPPFFGWSRFIPEGL", "QCSCGPDWYTVGTKYRSESYTWFLFIFCFIVPLSLICFSYTQLLRALKAVAAQQQESATT"]
    p3 = needleman_wunsch_profile(p1[0], p1[1], BLOSUM62)
    print(p3)
    #p4 = needleman_wunsch_profile(p3, p2[0], BLOSUM62)
    # print(p4)
    # p5 = needleman_wunsch_profile(p2[2], p2[1], BLOSUM62)
    # p6 = needleman_wunsch_profile(p4, p5, BLOSUM62)
    # print(p6)