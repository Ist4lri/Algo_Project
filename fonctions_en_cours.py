
from constants import *
from numpy import array, delete, zeros
from upgma import Upgma

class Alignement() :
    
    def __init__(self) :
        self.upgma = Upgma()
        self.clades_names = []
        self.matrice_distance = []
        self.temp_matrice = []
        self.temp_clades_names = []
        self.temp_tree = []
        self.freq_matrix = []
        self.dir_matrix = []


    def calculate_frequencies(self, seqs):
        """Fonction pour calculer les fréquences des acides aminés à chaque position

        Args:
            seqs (_type_): _description_
        """
        seq_len = len(seqs[0])
        self.freq_matrix = zeros((seq_len, 20))
        for seq in seqs:
            for i, aa in enumerate(seq):
                if aa == "-":
                    continue
                aa_index = "ARNDCQEGHILKMFPSTWY".index(aa) 
                self.freq_matrix[i][aa_index] += 1
        # Divise par la longueur :
        self.freq_matrix /= len(seqs)
        return self.freq_matrix


    def matrix_direction(self, seq1, seq2):
        self.dir_matrix = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
        # Remplissage de la matrice de direction
        for i in range(1, len(seq1)+1):
            for j in range(1, len(seq2)+1):
                if seq1[i-1] == seq2[j-1]:
                    self.dir_matrix[i][j] = 2  # Diagonale
                else:
                    max_val = max(self.dir_matrix[i-1][j], self.dir_matrix[i][j-1], self.dir_matrix[i-1][j-1])
                    if max_val == self.dir_matrix[i-1][j]:
                        self.dir_matrix[i][j] = 1  # Gauche
                    elif max_val == self.dir_matrix[i][j-1]:
                        self.dir_matrix[i][j] = 3  # Haut
                    else:
                        self.dir_matrix[i][j] = 2  # Diagonale

    def traceback(self, score_matrix, direction_matrix, seq1, seq2, gap="-"):
        """Fonction pour retracer l'alignement final à partir de la matrice de score et de la matrice de directions

        Args:
            score_matrix (_type_): _description_
            direction_matrix (_type_): _description_
            seq1 (_type_): _description_
            seq2 (_type_): _description_
            gap (str, optional): Char of the gap. Defaults to "-".

        Returns:
            _type_: _description_
        """
        self.freq_matrix[]
        i, j = len(seq1), len(seq2)
        aligned_seq1, aligned_seq2 = "", ""
        while i > 0 or j > 0:
            if direction_matrix[i][j] == "diag":
                aligned_seq1 = seq1[i-1] + aligned_seq1
                aligned_seq2 = seq2[j-1] + aligned_seq2
                i -= 1
                j -= 1
            elif direction_matrix[i][j] == "left":
                aligned_seq1 = gap + aligned_seq1
                aligned_seq2 = seq2[j-1] + aligned_seq2
                j -= 1
            else:
                aligned_seq1 = seq1[i-1] + aligned_seq1
                aligned_seq2 = gap + aligned_seq2
                i -= 1
        return aligned_seq1, aligned_seq2



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

