from constants import *
from numpy import array, delete, hstack, vstack

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
    counter_col = 0
    counter_lig = 0
    list_name_clades = []
    fasta_dict = open_fasta(filename)
    index_dict = len(fasta_dict)
    matrix_distance = array([[0] * (index_dict)
                            for _ in range(index_dict)])
    for i in range(len(matrix_distance)) :
        list_name_clades.append("seq"+str(i+1))
    for header in fasta_dict.keys():
        counter_col = 0
        seq_one = fasta_dict[header]
        for key, seq_two in fasta_dict.items():
            if header != key:
                matrix_distance[counter_lig, counter_col] = needleman_wunsch(
                    seq_one, seq_two, BLOSUM62)
            counter_col += 1
        counter_lig += 1
    return matrix_distance, list_name_clades


def find_min_upgma(matrix_distance):
    the_mini = {}
    mini = 100000
    the_mini[mini] = [0, 0]
    for x in range(0, len(matrix_distance)):
        for y in range(0, len(matrix_distance)):
            if x != y:
                temp_min = matrix_distance[x, y]
                if temp_min < mini:
                    del the_mini[mini]
                    mini = temp_min
                    the_mini[mini] = [x, y]
                    closer_clades = [x, y]
    return closer_clades

def paste_mini_upgma(matrix_dist, closer_clades, clades_names) :
    # on récupère les index des clades à regrouper trouvés plus haut
    clades_to_regroup = [clades_names[closer_clades[0]], clades_names[closer_clades[1]]]
    new_clades_names = []
    # les clades non concerné par la concaténation sont remis dans l'ordre dans la liste des noms de clades
    for seq_number in clades_names :
        if seq_number not in clades_to_regroup :
            new_clades_names.append(seq_number)
    # les clades regroupés sont ajoutés à la fin de la nouvelle liste (résultat = liste imbriquée des clades à associer)
    new_clades_names.append([clades_to_regroup])
    # on va crée une ligne/colonne qui aura les valeurs moyennes des deux lignes à supprmier
    new_col_of_matrix = []
    line_to_paste_1 = matrix_dist[closer_clades[0]]
    line_to_paste_2 = matrix_dist[closer_clades[1]]
    for i in range(len(line_to_paste_1)) :
        new_col_of_matrix.append((line_to_paste_1[i] + line_to_paste_2[i]) / 2)
    print(new_col_of_matrix)
    # Comme on va ajouter la colonne à la suite du tableau, il faut qu'on rajoute un zéro à la ligne (face à elle même)
    new_line_of_matrix = new_col_of_matrix
    new_line_of_matrix.append(0)
    print(new_line_of_matrix)
    # création de la nouvelle matrice :
    new_matrix_dist = matrix_dist
    print(new_matrix_dist)
    # on ajoute la colonne à la nouvelle matrice :
    # _______________________________
    # A PARTIR DE LÀ IL Y A ERREUR :
    new_matrix_dist = hstack([new_matrix_dist,new_col_of_matrix])
    # on ajoute la nouvelle ligne à la nouvelle matrice
    new_matrix_dist = vstack((new_matrix_dist, new_line_of_matrix))
    # on supprime les listes et colonnes correspondant au index des clades à regrouper
    new_matrix_dist = delete(new_matrix_dist, closer_clades[0], 0)
    new_matrix_dist = delete(new_matrix_dist, closer_clades[1], 0)
    new_matrix_dist = delete(new_matrix_dist, 0, closer_clades[0])
    new_matrix_dist = delete(new_matrix_dist, 0, closer_clades[1])
    print(new_matrix_dist)
    return new_matrix_dist, new_clades_names

def UPGMA(filename) :
    matrix_dist, clades_names = get_all_max_score(filename)
    closer_clades = find_min_upgma(matrix_dist)
    new_dist_matrix, new_clades_names = paste_mini_upgma(matrix_dist, closer_clades, clades_names)
    # mettre new_matrix avec nouveaux clades, si new_matrice de taille suffisante -> refaire upgma dessus, sinon arrêter
    return new_dist_matrix, new_clades_names

if __name__ == "__main__":
    print(UPGMA("opsines.fasta.txt"))
    
