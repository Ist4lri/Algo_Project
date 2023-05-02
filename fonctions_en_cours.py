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

    freq_matrix = np.zeros((21, seq_len))

    # Liste des acides aminés
    amino_acids = "CSTAGPDEQNHRKMILVWYF-"

    for i in range(num_seq):
        for j in range(seq_len) :
            indice_aa_matrice = amino_acids.find(seq_list[i][j])
            freq_matrix[indice_aa_matrice][j] += 1
    freq_matrix /= num_seq
    return freq_matrix


def needleman_wunsch_profile(seq1 : list, seq2 : list, matrix : list):
    """Cette fonction prend en entrée deux profils
    de séquences multiples, une matrice de score,
    et une valeur de gap (par défaut -4) et
    renvoie un nouveau profil aligné"""

    if type(seq1) is not list:
        seq1 = [seq1]
    if  type(seq2) is not list:
        seq2 = [seq2]

    # On trouve la séquence le profil le plus court (+ petite en ligne, + grande en colonnes)
    if len(seq1[0]) > len(seq2[0]):
        nb_col = len(seq1[0]) +1
        nb_ligne = len(seq2[0]) +1
        shorter_profile = seq2
        longer_profile = seq1
        freq_matrix1 = freq_matrix(seq1)
        freq_matrix2 = freq_matrix(seq2)
    else :
        nb_col = len(seq2[0]) +1
        nb_ligne = len(seq1[0]) +1
        shorter_profile = seq1
        longer_profile = seq2
        freq_matrix2 = freq_matrix(seq1)
        freq_matrix1 = freq_matrix(seq2)

    # définition du gap :
    gap = np.zeros(21)
    gap[-1] = 1
    # Initialisation de la matrice des scores et de direction
    score_matrix = np.zeros((nb_ligne, nb_col)) # m = nb de lignes, n = nb de colonnes
    direction_matrix = np.zeros((nb_ligne, nb_col))
    
    # On rempli la première colonne en fonction des gaps:
    for i in range(1, nb_ligne):
        # on prend le score de la case de gauche,
        # on l'additionne à la np.array([0,0,0, ...,1]) @ BLOSUM62 @ fréq
        score_matrix[i][0] = score_matrix[i-1][0] + (gap @ matrix @ freq_matrix2)[i-1]
        direction_matrix[i][0] = 1
    # On rempli la prmière ligne en fonction des gaps:
    for j in range(1, nb_col):
        # on prend le score de la case du haut,
        # on l'additionne à la fréquence d'un gap * BLOSUM62 * la fréquence observée sur la position
        score_matrix[0][j] = score_matrix[0][j-1] + (gap @ matrix @ freq_matrix1)[j-1]
        direction_matrix[0][j] = 2

    print("\nscore_matrix ap : \n", score_matrix)
    
    # On change de sens pour pouvoir faire le calcul matriciel
    freq_matrix1 = np.transpose(freq_matrix1)
    freq_matrix2 = np.transpose(freq_matrix2)
    

    # Remplissage de la matrice des scores et de direction
    for i in range(1, nb_col): # taille de la seq la plus longue (seq1)
        for j in range(1, nb_ligne): # taille de la seq la plus courte (seq2)
            # Calcul du score pour un alignement diagonal
            diag_score = score_matrix[j-1][i-1] + (freq_matrix1[i-1] @ matrix @ freq_matrix2[j-1])
            
            # Calcul du score pour un alignement horizontal
            left_score = score_matrix[j][i-1] + (freq_matrix2[j-1] @ matrix @ gap)
            #left_score = score_matrix[i-1][j] + (freq_matrix2[j-1] @ matrix @ gap)
            
            # Calcul du score pour un alignement vertical
            up_score = score_matrix[j-1][i] + (freq_matrix1[i-1] @ matrix @ gap)
            #up_score = score_matrix[i][j-1] + (freq_matrix1[i-1] @ matrix @ gap)
            
            score_matrix[j][i], direction_matrix[j][i] = max((diag_score, 0), (up_score, 1), (left_score, 2))

        # aller à gauche met un gap sur la seq de gauche, aller en haut met un gap sur la seq du haut.
    print("score_matrice après calcul : \n", score_matrix)
    print("distance_matrice après calcul : \n", direction_matrix)
    # Construction du nouvel alignement
    aligned_seq1 = []
    aligned_seq2 = []

    # Pour chaque séquence de la liste/du profile
    for k in longer_profile:
        i = nb_ligne - 1
        j = nb_col - 1
        #print('\n new_aa_to_join: ',new_aa_to_join)
        new_aa_to_join = ""
        # On va chercher dans la matrice de direction
        while i > 0 or j > 0:
            direction = direction_matrix[i][j]
            # Si le score vient de la diagonale :
            if direction == 0:
                new_letter = k[j-1]
                new_aa_to_join += new_letter
                i -= 1
                j -= 1
            # Si le score vient d'en haut
            elif direction == 1:
                new_aa_to_join += "-"
                i -= 1
            # Si le score vient de la gauche
            elif direction == 2:
                new_letter = k[j-1]
                new_aa_to_join += new_letter
                j -= 1
        aligned_seq1.append(new_aa_to_join[::-1])
        print("aligned_seq1 : ", aligned_seq1)

    # Pour chaque séquence dans le profile
    for k in shorter_profile:
        i = nb_ligne - 1
        j = nb_col - 1
        new_aa_to_join2 = ""
        # Tant qu'il reste des AA ou nt
        while i > 0 or j > 0:
            # On regarde dans la matrice de direction
            direction = direction_matrix[i][j]
            # Si le score vient de la diagonale
            if direction == 0:
                new_letter = k[i-1]
                new_aa_to_join2 += new_letter
                i -= 1
                j -= 1
            # Si le score vient d'en haut
            elif direction == 1:
                new_letter = k[i-1]
                new_aa_to_join2 += new_letter
                i -= 1
            # Si le score vient de la gauche
            elif direction == 2:
                new_aa_to_join2 += "-"
                j -= 1
        aligned_seq2.append(new_aa_to_join2[::-1])

    new_alignment = []

    for elem in aligned_seq1:
        new_alignment.append(elem)
    for elem in aligned_seq2:
        new_alignment.append(elem)

    print("\n NEW ALIGNEMENTS: ", new_alignment, "\n __________________________________________________________ \n")

    return new_alignment


if __name__ == "__main__":

#     ptest_2sequence = ["MSKMSEEEEFLLFKNISLVGPWDGPQYHLAPVWAFHLQAVFMGFVFFVGTPLNATVLVAT\
# LRYRKLRQPLNYILVNVSLGGFIYCIFSVFIVFITSCYGYFVFGRHVCALEAFLGCTAGL\
# VTGWSLAFLAFERYIIICKPFGNFRFSSKHALMVVVATWTIGIGVSIPPFFGWSRFVPEG\
# LQCSCGPDWYTVGTKYYSEYYTWFLFIFCYIVPLSLICFSYSQLLGALRAVAAQQQESAS\
# TQKAEREVSHMVVVMVGSFCLCYTPYAALAMYIVNNRNHGVDLRLVTIPAFFSKSACVYN\
# PIIYCFMNKQFRACIMEMVCGKPMTDESELSSSQKTEVSTVSSSQVGPN","MDAWAVQFGNASKVSPFEGEQYHIAPKWAFYLQAAFMGFVFIVGTPMNGIVLFVTMKYKK\
# LRQPLNYILVNISLAGFIFDTFSVSQVFVCAARGYYFLGYTLCAMEAAMGSIAGLVTGWS\
# LAVLAFERYVVICKPFGSFKFGQGQAVGAVVFTWIIGTACATPPFFGWSRYIPEGLGTAC\
# GPDWYTKSEEYNSESYTYFLLITCFMMPMTIIIFSYSQLLGALRAVAAQQAESESTQKAE\
# REVSRMVVVMVGSFVLCYAPYAVTAMYFANSDEPNKDYRLVAIPAFFSKSSCVYNPLIYA"]

#     p_test = needleman_wunsch_profile(ptest_2sequence[0], ptest_2sequence[1], BLOSUM62)
#     print(p_test)
#     #p4 = needleman_wunsch_profile(p3, p2[0], BLOSUM62)
#     # print(p4)
#     # p5 = needleman_wunsch_profile(p2[2], p2[1], BLOSUM62)
#     # p6 = needleman_wunsch_profile(p4, p5, BLOSUM62)
#     # print(p6)

    bovin = ["MSKMSEEEEFLLFKNISLVGPWDGPQYHLAPVWAFHLQAVFMGFVFFVGTPLNATVLVATLRYRKLRQPLNYILVNVSLGGFIYCIFSVFIVFITSCYGYFVFGRHVCALEAFLGCTAGLVTGWSLAFLAFERYIIICKPFGNFRFSSKHALMVVVATWTIGIGVSIPPFFGWSRFVPEGLQCSCGPDWYTVGTKYYSEYYTWFLFIFCYIVPLSLICFSYSQLLGALRAVAAQQQESASTQKAEREVSHMVVVMVGSFCLCYTPYAALAMYIVNNRNHGVDLRLVTIPAFFSKSACVYNPIIYCFMNKQFRACIMEMVCGKPMTDESELSSSQKTEVSTVSSSQVGPN"]
    
    Homo = ["MRKMSEEEFYLFKNISSVGPWDGPQYHIAPVWAFYLQAAFMGTVFLIGFPLNAMVLVATL\
    RYKKLRQPLNYILVNVSFGGFLLCIFSVFPVFVASCNGYFVFGRHVCALEGFLGTVAGLV\
    TGWSLAFLAFERYIVICKPFGNFRFSSKHALTVVLATWTIGIGVSIPPFFGWSRFIPEGL\
    QCSCGPDWYTVGTKYRSESYTWFLFIFCFIVPLSLICFSYTQLLRALKAVAAQQQESATT\
    QKAEREVSRMVVVMVGSFCVCYVPYAAFAMYMVNNRNHGLDLRLVTIPSFFSKSACIYNP\
    IIYCFMNKQFQACIMKMVCGKAMTDESDTCSSQKTEVSTVSSTQVGPN"]

    Macaca = ["MRKMSEEEEFYLFKNLSSVKPWDGPQYHIAPVWAFYLQAAFMGTVFLAGFPLNAMVLVAT\
    VRYKKLRQPLNYILVNVSFGGFLLCIFSVFPVFVNSCKGYFVFGRHVCAFEAFLGTVAGL\
    VTGWSLAFLAFERYIVICKPFGNFRFSSKHALTVVLATWTIGIGVSIPPFFGWSRFIPEG\
    LQCSCGPDWYTVGTKYRSESYTWFLFIFCFIVPLSLICFSYTQLLRALKAVAAQQQESAT\
    TQKAEREVSRMVVVMVGSFCVCYVPYAAFAMYMVNNRNHGLDLRLVTIPAFFSKSACIYN\
    PIIYCFMNKQFQAHIMKMVCGKAMTDESDISSSQKTEVSTVSSSQVGPN"]

    Macropus = ["MSGDEEFYLFKNISSVGPWDGPQYHIAPAWAFHCQTVFMGFVFFAGTPLNAVVLIATFRYKKLRQPLNYILVNISLAGFIYCIFSVFTVFISSSQGYFIFGRHVCAMEGFLGSVAGLVTGWSLAFLAFERFIVICKPFGNFRFNSKHSMMVVLATWVIGIGVSIPPFFGWSRYIPEGLQCSCGPDWYTVGTKYRSEYYTWFLFILCFIMPLSLICFSYSQLLGALRAVAAQQQESATTQKAEREVSRMVVMMVGSFCLCYVPYAALAMYMVNNRNHGIDLRLVTIPAFFSKSSCVYNPIIYCFMNKQFHACIMEMVCRKPMTDDSEASSSQKTEVSTVSSSQVGPS"]
    
    Profile_test1 = [bovin, Macropus]
    #Profile_test2 = [Macaca, Macropus]
    
    p_test1 = needleman_wunsch_profile(Profile_test1[0], Profile_test1[1], BLOSUM62)
    print(p_test1)
    # p_test2 = needleman_wunsch_profile(Profile_test2[0], Profile_test2[1], BLOSUM62)
    # p_test3 = needleman_wunsch_profile(p_test1, p_test2, BLOSUM62)
    # print(p_test3)
