
from constants import *
from numpy import array, delete, zeros
from upgma import Upgma

# class Alignement() :
    
#     def __init__(self) :
#         self.upgma = Upgma()
#         self.clades_names = []
#         self.matrice_distance = []
#         self.temp_matrice = []
#         self.temp_clades_names = []
#         self.temp_tree = []
#         self.freq_matrix = []
#         self.dir_matrix = []


#     def calculate_frequencies(self, seqs):
#         """Fonction pour calculer les fréquences des acides aminés à chaque position

#         Args:
#             seqs (_type_): _description_
#         """
#         seq_len = len(seqs[0])
#         self.freq_matrix = zeros((seq_len, 20))
#         for seq in seqs:
#             for i, aa in enumerate(seq):
#                 if aa == "-":
#                     continue
#                 aa_index = "ARNDCQEGHILKMFPSTWY".index(aa) 
#                 self.freq_matrix[i][aa_index] += 1
#         # Divise par la longueur :
#         self.freq_matrix /= len(seqs)
#         return self.freq_matrix


#     def matrix_direction(self, seq1, seq2):
#         self.dir_matrix = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
#         # Remplissage de la matrice de direction
#         for i in range(1, len(seq1)+1):
#             for j in range(1, len(seq2)+1):
#                 if seq1[i-1] == seq2[j-1]:
#                     self.dir_matrix[i][j] = 2  # Diagonale
#                 else:
#                     max_val = max(self.dir_matrix[i-1][j], self.dir_matrix[i][j-1], self.dir_matrix[i-1][j-1])
#                     if max_val == self.dir_matrix[i-1][j]:
#                         self.dir_matrix[i][j] = 1  # Gauche
#                     elif max_val == self.dir_matrix[i][j-1]:
#                         self.dir_matrix[i][j] = 3  # Haut
#                     else:
#                         self.dir_matrix[i][j] = 2  # Diagonale

#     # def traceback(self, score_matrix, direction_matrix, seq1, seq2, gap="-"):
#     #     """Fonction pour retracer l'alignement final à partir de la matrice de score et de la matrice de directions

#     #     Args:
#     #         score_matrix (_type_): _description_
    #         direction_matrix (_type_): _description_
    #         seq1 (_type_): _description_
    #         seq2 (_type_): _description_
    #         gap (str, optional): Char of the gap. Defaults to "-".

    #     Returns:
    #         _type_: _description_
    #     """
    #     self.freq_matrix[]
    #     i, j = len(seq1), len(seq2)
    #     aligned_seq1, aligned_seq2 = "", ""
    #     while i > 0 or j > 0:
    #         if direction_matrix[i][j] == "diag":
    #             aligned_seq1 = seq1[i-1] + aligned_seq1
    #             aligned_seq2 = seq2[j-1] + aligned_seq2
    #             i -= 1
    #             j -= 1
    #         elif direction_matrix[i][j] == "left":
    #             aligned_seq1 = gap + aligned_seq1
    #             aligned_seq2 = seq2[j-1] + aligned_seq2
    #             j -= 1
    #         else:
    #             aligned_seq1 = seq1[i-1] + aligned_seq1
    #             aligned_seq2 = gap + aligned_seq2
    #             i -= 1
    #     return aligned_seq1, aligned_seq2



# # Fonction pour calculer le score de substitution entre deux profils de fréquences
# def score_frequency_matrix(f1, f2, gap_penalty=-5, matrix="BLOSUM62"):
#     score_matrix = np.zeros((f1.shape[0]+1, f2.shape[0]+1))
#     for i in range(1, f1.shape[0]+1):
#         for j in range(1, f2.shape[0]+1):
#             match_score = np.sum(f1[i-1]*f2[j-1]*BLOSUM62)
#             gap_penalty_1 = score_matrix[i-1][j] + gap_penalty
#             gap_penalty_2 = score_matrix[i][j-1] + gap_penalty
#             match_mismatch_score = score_matrix[i-1][j-1] + match_score
#             score_matrix[i][j] = max(gap_penalty_1, gap_penalty_2, match_mismatch_score)
#     return score_matrix

# # Fonction pour retracer l'alignement final à partir de la matrice de score et de la matrice de directions
# def traceback(score_matrix, direction_matrix, seq1, seq2, gap="-"):
#     i, j = len(seq1), len(seq2)
#     aligned_seq1, aligned_seq2 = "", ""
#     while i > 0 or j > 0:
#         if direction_matrix[i][j] == "diag":
#             aligned


# # Fonction pour retracer l'alignement final à partir de la matrice de score et de la matrice de directions
# def traceback(score_matrix, direction_matrix, seq1, seq2, gap="-"):
#     i, j = len(seq1), len(seq2)
#     aligned_seq1, aligned_seq2 = "", ""
#     while i > 0 or j > 0:
#         if direction_matrix[i][j] == "diag":
#             aligned_seq1 = seq1[i-1] + aligned_seq1
#             aligned_seq2 = seq2[j-1] + aligned_seq2
#             i -= 1
#             j -= 1
#         elif direction_matrix[i][j] == "left":
#             aligned_seq1 = gap + aligned_seq1
#             aligned_seq2 = seq2[j-1] + aligned_seq2
#             j -= 1
#         else:
#             aligned_seq1 = seq1[i-1] + aligned_seq1
#             aligned_seq2 = gap + aligned_seq2
#             i -= 1
#     return aligned_seq1, aligned_seq2

# # Fonction récursive pour faire l'alignement multiple en suivant l'arbre guide
# def align_sequences(node, alignments_dict):
#     if "clades" in node:
#         clade1, clade2 = node["clades"]
#         align_sequences(alignments_dict[clade1], alignments_dict)
#         align_sequences(alignments_dict[clade2], alignments_dict)
#         seqs1, seqs2 = alignments_dict[clade1], alignments_dict[clade2]
#         freq_matrix1 = calculate_frequencies(seqs1)
#         freq_matrix2 = calculate_frequencies(seqs2)
#         score_matrix = score_frequency_matrix(freq_matrix1, freq_matrix2)
#         direction_matrix = np.zeros_like(score_matrix, dtype=np.object)
#         for

import numpy as np

# def freq_matrix(seq_list):
#     """Cette fonction prend une liste de séquences d'acides aminés et retourne une matrice de fréquence"""
#     num_seq = len(seq_list)
#     seq_len = len(seq_list[0])
#     freq_matrix = np.zeros((seq_len, 20))
#     for i in range(seq_len):
#         for j in range(num_seq):
#             if seq_list[j][i] != "-":
#                 freq_matrix[i][ord(seq_list[j][i])-65] += 1
#     freq_matrix /= num_seq
#     return freq_matrix

def freq_matrix(seq_list):
    """Cette fonction prend une liste de séquences d'acides aminés et retourne une matrice de fréquence"""
    
    # Liste des acides aminés
    amino_acids = "CSTAGPDEQNHRKMILVWYF-"
    
    # Initialisation de la matrice avec des zéros
    num_seq = len(seq_list)
    seq_len = len(seq_list[0])
    freq_matrix = np.zeros((seq_len, len(amino_acids)))
    
    # Remplissage de la matrice de fréquences
    for i in range(seq_len):
        for j in range(num_seq):
            indice_aa_dans_matrice = amino_acids.find(seq_list[j][i])
            freq_matrix[i][indice_aa_dans_matrice] += 1
    freq_matrix /= num_seq
    return freq_matrix


def needleman_wunsch_profile(seq1 : list, seq2 : list, matrix : list, gap=-4):
    """Cette fonction prend en entrée deux profils de séquences multiples, une matrice de score, et une valeur de gap (par défaut -4) et renvoie un nouveau profil aligné"""

    if type(seq1) is not list :
        seq1 = [seq1]
    if type(seq2) is not list :
        seq2 = [seq2]
    len_seq1 = len(seq1[0]) + 1
    len_seq2 = len(seq2[0]) + 1
    
    freq_matrix1 = freq_matrix(seq1)
    freq_matrix2 = freq_matrix(seq2)
    
    # Initialisation de la matrice des scores et de direction
    score_matrix = np.zeros((len_seq1, len_seq2))
    direction_matrix = np.zeros((len_seq1, len_seq2))
    for i in range(1, len_seq1):
        score_matrix[i][0] = score_matrix[i-1][0] + gap
        direction_matrix[i][0] = 1
    for j in range(1, len_seq2):
        score_matrix[0][j] = score_matrix[0][j-1] + gap
        direction_matrix[0][j] = 2
    
    # Remplissage de la matrice des scores et de direction
    for i in range(1, len_seq1):
        for j in range(1, len_seq2):
            # Calcul du score pour un alignement diagonal
            diag_score = score_matrix[i-1][j-1] + np.sum(freq_matrix1[i-1] * matrix * freq_matrix2[j-1])
            # Calcul du score pour un alignement vertical
            up_score = score_matrix[i-1][j] + gap
            # Calcul du score pour un alignement horizontal
            left_score = score_matrix[i][j-1] + gap
            # Choix du score maximum et de la direction correspondante
            score_matrix[i][j], direction_matrix[i][j] = max((diag_score, 0), (up_score, 1), (left_score, 2))
    
    print("\nMatrice des scores : \n", score_matrix)
    #print("\nMatrice des directions : \n", direction_matrix, "\n")
    
    # Construction du nouvel alignement
    aligned_seq1 = []
    aligned_seq2 = []
    
    for k in seq1 :
        new_seq_in_progress1 = ""
        new_seq_in_progress2 = ""
        i = len_seq1 - 1
        j = len_seq2 - 1
        while i > 0 or j > 0:
            direction = direction_matrix[i][j]
            new_letter = k[i-1]
            if direction == 0:
                new_seq_in_progress1 += new_letter
                i -= 1
                j -= 1
            elif direction == 1:
                new_seq_in_progress1 += new_letter
                i -= 1
            else:
                new_seq_in_progress1 += "-"
                j -= 1
        aligned_seq1.append(new_seq_in_progress1[::-1])

    for k in seq2 :
        new_seq_in_progress1 = ""
        new_seq_in_progress2 = ""
        i = len_seq1 - 1
        j = len_seq2 - 1
        while i > 0 or j > 0:
            direction = direction_matrix[i][j]
            new_letter = k[i-1]
            if direction == 0:
                new_seq_in_progress2 += new_letter
                i -= 1
                j -= 1
            elif direction == 1:
                new_seq_in_progress2 += "-"
                i -= 1
            else:
                new_seq_in_progress2 += new_letter
                j -= 1
        aligned_seq2.append(new_seq_in_progress2[::-1])
    
    print("\n aligned_seq 1 : ", aligned_seq1)
    print("\n aligned_seq2 : ", aligned_seq2)
    
    # Convertir les listes de séquences alignées en une liste de séquences alignées unique
    new_alignment = []
    for i in aligned_seq1:
        new_alignment.append(i)
    for i in aligned_seq2:
        new_alignment.append(i)

    return new_alignment


if __name__ == "__main__":
    p1 = ["MSKMSEEEEFLLFKNISLVGPWDGPQYHLAPVWAFHLQAVFMGFV", "MRKMSEEEFYLFKNISSVGPWDGPQYHIAPVWAFYLQAAFMGTVFLIGFPLNAMVLVATL"]
    p2 = ["MSKMSEEEEFLLFKNISLVGPWDGPQYHLAPVWAFHLQAVFMGFV", "MRKMSEEEFYLFKNISSVGPWDGPQYHIAPVWAFYLQAAFMGTVFLIGFPLNAMVLVATL"]
    # p3 = ["MRKMSEEEFYLFKNISSVGPWDGPQYHIAPVWAFYLQAAFMGTVFLIGFPLNAMVLVATL", "RYKKLRQPLNYILVNVSFGGFLLCIFSVFPVFVASCNGYFVFGRHVCALEGFLGTVAGLV",
    #       "TGWSLAFLAFERYIVICKPFGNFRFSSKHALTVVLATWTIGIGVSIPPFFGWSRFIPEGL", "QCSCGPDWYTVGTKYRSESYTWFLFIFCFIVPLSLICFSYTQLLRALKAVAAQQQESATT"]
    p3 = needleman_wunsch_profile(p1[0], p2[0], BLOSUM62)
    p4 = needleman_wunsch_profile(p1[1], p2[1], BLOSUM62)
    p5 = needleman_wunsch_profile(p3, p4, BLOSUM62)
    print(p5[0])
    print(p5[1])
    print(p5[2])
    print(p5[3])
    # p4 = needleman_wunsch_profile(p3, p2[0], BLOSUM62)
    # p5 = needleman_wunsch_profile(p2[2], p2[1], BLOSUM62)
    # p6 = needleman_wunsch_profile(p4, p5, BLOSUM62)
    # print(p6)