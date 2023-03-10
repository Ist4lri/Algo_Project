from constants import *
from numpy import array, delete, zeros
from upgma import Upgma

class Alignement() :
    
    def __init__(self) :
        self.upgma = Upgma()
        self.clades_names = [] #Liste des noms de séquences (en seqX)
        self.matrice_distance = [] #Matrice de distance
        self.freq_matrix = [] #Matrice de fréquences
        self.dir_matrix = [] #Matrice des directions (Après 2e NW)
        self.dict_tree = {} #Dico après premier UPGMA arbre {'seqX' : distance}
        self.dict_sequence_tree = {} #Dico avec les séquences pour alignement {'seqX' :  AGTCGAT...}
        self.group_sequence_align = [] #Multiple alignement
        self.second_groupe_align = [] #Si y'a un autre multiple alignement
        
    
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
        # To UPGMA :
        self.to_upgma(self.matrice_distance, self.clades_names)
        
        
    def to_upgma(self, distance_matrix, clades_names):
        self.dict_tree  = self.upgma.tree_with_upgma(distance_matrix, clades_names)
        self.multiple_alignement()

    def freq_matrix(self, seq_list):
        """Cette fonction prend une liste de séquences d'acides aminés et retourne une matrice de fréquence"""
        num_seq = len(seq_list)
        seq_len = len(seq_list[0])
        freq_matrix = np.zeros((seq_len, 20))
        for i in range(seq_len):
            for j in range(num_seq):
                if seq_list[j][i] != "-":
                    freq_matrix[i][ord(seq_list[j][i])-65] += 1
        freq_matrix /= num_seq
        return freq_matrix

    def needleman_wunsch_profile(self, seq1, seq2, matrix, gap=-4):
        """Cette fonction prend en entrée deux profils de séquences multiples, une matrice de score, et une valeur de gap (par défaut -4) et renvoie un nouveau profil aligné"""
        
        # Calcul des matrices de fréquence pour les deux profils
        freq_matrix1 = self.freq_matrix(seq1)
        freq_matrix2 = self.freq_matrix(seq2)
        
        # Initialisation de la matrice des scores et de direction
        m = len(seq1[0]) + 1
        n = len(seq2[0]) + 1
        score_matrix = zeros((m, n))
        direction_matrix = zeros((m, n))
        for i in range(1, m):
            score_matrix[i][0] = score_matrix[i-1][0] + gap
            direction_matrix[i][0] = 1
        for j in range(1, n):
            score_matrix[0][j] = score_matrix[0][j-1] + gap
            direction_matrix[0][j] = 2
        
        # Remplissage de la matrice des scores et de direction
        for i in range(1, m):
            for j in range(1, n):
                # Calcul du score pour un alignement diagonal
                diag_score = score_matrix[i-1][j-1] + sum(freq_matrix1[i-1] * matrix * freq_matrix2[j-1])
                # Calcul du score pour un alignement vertical
                up_score = score_matrix[i-1][j] + gap
                # Calcul du score pour un alignement horizontal
                left_score = score_matrix[i][j-1] + gap
                # Choix du score maximum et de la direction correspondante
                score_matrix[i][j], direction_matrix[i][j] = max((diag_score, 0), (up_score, 1), (left_score, 2))
        
        # Construction du nouvel alignement
        aligned_seq1 = []
        aligned_seq2 = []
        i = m - 1
        j = n - 1
        while i > 0 or j > 0:
            direction = direction_matrix[i][j]
            if direction == 0:
                aligned_seq1.append("".join(seq1[k][i-1] for k in range(len(seq1))))
                aligned_seq2.append("".join(seq2[k][j-1] for k in range(len(seq2))))
                i -= 1
                j -= 1
            elif direction == 1:
                aligned_seq1.append("".join(seq1[k][i-1] for k in range(len(seq1))))
                aligned_seq2.append("-" * len(seq2[0]))
                i -= 1
            else:
                aligned_seq1.append("-" * len(seq1[0]))
                aligned_seq2.append("".join(seq2[k][j-1] for k in range(len(seq2))))
                j -= 1
        
        # Inverse les séquences alignées
        aligned_seq1.reverse()
        aligned_seq2.reverse()
        
        # Convertir les listes de séquences alignées en une liste de séquences alignées unique
        num_seq = len(seq1) + len(seq2)
        new_alignment = []
        for i in range(num_seq):
            new_alignment.append(aligned_seq1[i] + aligned_seq2[i-len(seq1)])
        
        return new_alignment
    
    def multiple_alignement(self):
        list_seq = []
        list_dist = []
        for k, v in sorted(self.dict_tree.items(), key=lambda x: x[1]):
            del self.dict_tree[k]
            self.dict_tree[k] = v
        for key, value in self.dict_tree.items():
            list_seq.append(key)
            list_dist.append(value)
        
        for i in range(len(list_dist)):
            if i+1 < len(list_dist):
                if list_dist[i] == list_dist[i+1]:
                    self.needleman_wunsch_profile(list_dist)
                    print("A comparer dans []*, puis comparer à []")
                else:
                    print("(A vérifier si la seq a déjà été utilisée...)A comparer avec le []")
            else:
                print("ok")
            
            if i+1 != i+2:
                print("bite")
                #Merge les seq i et i+1
                #Calcul les freq entre seq i et i+1 de cette matrice
                #alignement needleman, matrice direction, qui ressort tuple ou liste de seqs (même taille avec gap)
                #update dico pour fusionner les 2
            else:
                print("Oui")
                #Merge i+1 et i+2
                #
        
            
    
 
if __name__ == "__main__":
    align = Alignement()
    print(align.calculate_frequencies(["ARNDCQEGHILKMFPSTWY", "ARNDCQEGHILKMFPSTWY", "ARNDCQEGHILKMFPSTWY"]))

