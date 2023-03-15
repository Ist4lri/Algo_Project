from constants import *
from numpy import array, delete, zeros, ones
from numpy import max as np_max
from upgma import Upgma


class Alignement():

    def __init__(self):
        self.upgma = Upgma()
        self.clades_names = []  # Liste des noms de séquences (en seqX)
        self.matrice_distance = []  # Matrice de distance
        self.dir_matrix = []  # Matrice des directions (Après 2e NW)
        self.dict_tree = {}  # Dico après premier UPGMA arbre {'seqX' : distance}
        # Dico avec les séquences pour alignement {'seqX' :  AGTCGAT...}
        self.dict_sequence_tree = {}
        self.alignement1 = []  # Premier alignement
        self.alignement2 = []  # Deuxième alignement si nécessaire
        self.freq_mat1 = []  # matrice de fréquences pour l'alignement1
        self.freq_mat2 = []  # matrice de fréquences pour l'alignement2

    def needleman_wunsch(self, seq1, seq2, matrix, gap=-1) -> int:
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

    def get_all_max_score(self, fasta_dict):
        """get the max score of all sequence

        Args:
            filename (str): path of the file

        Returns:
            list: array of array with all score
        """
        counter_col = 0
        counter_lig = 0
        index_dict = len(fasta_dict)
        # Put only zeros in a mat of the good size
        self.matrice_distance = array([[0] * (index_dict)
                                       for _ in range(index_dict+1)])
        # shorter names :
        for i in range(len(self.matrice_distance)+1):
            self.clades_names.append("seq"+str(i+1))
        # For each pair of seq :
        for header in fasta_dict.keys():
            counter_col = 0
            counter_lig += 1
            seq_one = fasta_dict[header]
            for key, seq_two in fasta_dict.items():
                if header != key:
                    # find the score between both seq and put it in the distance mat :
                    self.matrice_distance[counter_lig, counter_col] = self.needleman_wunsch(
                        seq_one, seq_two, BLOSUM62)
                counter_col += 1
        # To delete the first line of zeros :
        self.matrice_distance = delete(self.matrice_distance, 0, axis=0)
        # Transform Needleman-Wunsch matrix to an Upgma matrix :
        len_mat_dist = len(self.matrice_distance)
        print("=============AVANT MODIF=============")
        print(self.matrice_distance)
        self.matrice_distance = self.matrice_distance - \
            np_max(self.matrice_distance)*ones((len_mat_dist, len_mat_dist))
        for line in range(len(self.matrice_distance)):
            for col in range(len(self.matrice_distance)):
                if self.matrice_distance[line, col] != 0:
                    self.matrice_distance[line, col] = - \
                        self.matrice_distance[line, col]
        print("=============APRES MODIF=============")
        print(self.matrice_distance)
        #self.matrice_distance = self.upgma.transform_mat_dist(self.matrice_distance)
        # To UPGMA :
        self.upgma_to_multiple_align(self.matrice_distance, self.clades_names)

    def upgma_to_multiple_align(self, distance_matrix, clades_names):
        self.dict_tree = self.upgma.tree_with_upgma(
            distance_matrix, clades_names)
        self.multiple_alignement()

    def freq_matrix_calc(self, seq_list):
        """Cette fonction prend une liste de séquences d'acides aminés et retourne une matrice de fréquence"""
        # Liste des acides aminés
        amino_acids = "CSTAGPDEQNHRKMILVWYF-"
        # Initialisation de la matrice avec des zéros
        num_seq = len(seq_list)
        seq_len = len(seq_list[0])
        freq_matrix = zeros((seq_len, len(amino_acids)))
        # Remplissage de la matrice de fréquences
        for i in range(seq_len):
            for j in range(num_seq):
                indice_aa_dans_matrice = amino_acids.find(seq_list[j][i])
                freq_matrix[i][indice_aa_dans_matrice] += 1
        freq_matrix /= num_seq
        return freq_matrix

    def needleman_wunsch_profile(self, seq1: list, seq2: list, matrix: list, gap=-4):
        """Cette fonction prend en entrée deux profils de séquences multiples, une matrice de score,
        et une valeur de gap (par défaut -4) et renvoie un nouveau profil aligné"""
        if type(seq1) is not list:
            seq1 = [seq1]
        if type(seq2) is not list:
            seq2 = [seq2]
        len_seq1 = len(seq1[0]) + 1
        len_seq2 = len(seq2[0]) + 1

        freq_mat1 = self.freq_matrix_calc(seq1)
        freq_mat2 = self.freq_matrix_calc(seq2)

        print("\nF_mat1 : ", freq_mat1, "\n")
        print("\nF_mat2 : ", freq_mat2, "\n")

        # Initialisation de la matrice des scores et de direction
        score_matrix = zeros((len_seq1, len_seq2))
        direction_matrix = zeros((len_seq1, len_seq2))
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
                print("\nA : ", seq1.index(seq1[i-1]))
                print("\nB : ", seq2.index(seq2[j-1]))
                print("\nC : ", matrix[seq1.index(
                    seq1[i-1]), seq2.index(seq2[j-1])])
                diag_score = sum(score_matrix[i-1][j-1] + freq_mat1[i-1] * matrix[seq1.index(
                    seq1[i-1]), seq2.index(seq2[j-1])] * freq_mat2[j-1])
                print("\nDiag_score : ", diag_score)
                # Calcul du score pour un alignement vertical
                up_score = score_matrix[i-1][j] + gap
                print("\nUp_score : ", up_score)
                # Calcul du score pour un alignement horizontal
                left_score = score_matrix[i][j-1] + gap
                print("\nLeft_score : ", left_score)
                # Choix du score maximum et de la direction correspondante
                score_matrix[i][j], direction_matrix[i][j] = max(
                    (diag_score, 0), (up_score, 1), (left_score, 2))

        # Construction du nouvel alignement
        aligned_seq1 = []
        aligned_seq2 = []

        for k in seq1:
            new_seq_in_progress1 = ""
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

        for k in seq2:
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

        # Convertir les listes de séquences alignées en une liste de séquences alignées unique
        new_alignment = []
        for i in aligned_seq1:
            new_alignment.append(i)
        for i in aligned_seq2:
            new_alignment.append(i)

        return new_alignment

    def multiple_alignement(self):
        """_summary_
        """
        list_seq = []
        list_dist = []
        # On trie le dico avec en clé les séquences et en valeur la distance (les plus proches en premiers)
        for k, v in sorted(self.dict_tree.items(), key=lambda x: x[1]):
            del self.dict_tree[k]
            self.dict_tree[k] = v
        # On place les clés et valeurs dans deux listes séparées (c'est trié) :
        for key, value in self.dict_tree.items():
            list_seq.append(key)
            list_dist.append(value)
        # On sort les séquences à merge du bon dico qui les contient :
        print(self.dict_tree)

        for i in range(len(list_dist)-1):
            print("\n Etape 0 : ", i)
            # Si les deux première distances (triées) sont identiques :
            if list_dist[i] == list_dist[i+1]:
                espece1 = self.dict_sequence_tree[list_seq[i]]
                espece2 = self.dict_sequence_tree[list_seq[i+1]]
                self.alignement1 = self.needleman_wunsch_profile(
                    espece1, espece2, BLOSUM62)
                i += 2
                print("\n Etape 1 : ", i)
                print(self.alignement1)
                if list_dist[i] == list_dist[i+1]:
                    espece3 = self.dict_sequence_tree[list_seq[i]]
                    espece4 = self.dict_sequence_tree[list_seq[i+1]]
                    i += 2
                    print("\n Etape 2 : ", i)
                    self.alignement2 = self.needleman_wunsch_profile(
                        espece3, espece4, BLOSUM62)
                if self.alignement2 != []:
                    self.alignement1 = self.needleman_wunsch_profile(
                        self.alignement1, self.alignement2, BLOSUM62)
                    self.alignement2 = []
            else:
                print("\n Etape 3 : ", i)
                espece_solo = self.dict_sequence_tree[list_seq[i]]
                self.alignement1
                print("bite")
        print("apres")
        return self.alignement1


if __name__ == "__main__":
    align = Alignement()
    print(align.multiple_alignement(
        ["ARNDCQEGHILKMFPSTWY", "ARNDCQEGHILKMFPSTWY", "ARNDCQEGHILKMFPSTWY"]))
