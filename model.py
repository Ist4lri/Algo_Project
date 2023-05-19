from arbre import *
from constants import *
from alignement import *
import argparse
from os import path
from sys import exit as sys_exit


def is_fasta(fastafile) : # takes in argument the name of the file
    """This function control if the file starts with a '>' used in fasta files headers."""
    if fastafile.split(".")[-1] in ["fasta","faa","fa"] : # extension is OK
        with open(fastafile, "r", encoding='UTF-8') as file :
            first_line = file.read(1) # take the first line of the file
            if first_line.split()[0] == ">" : # first char of first line
                return True # this looks like a fasta file
    return False


class Model():
    """Model Class"""

    def __init__(self):
        """Init of the class"""
        self.dict_of_seq = {}
        self.dict_real_name = {}
        self.alignement = Alignement()

    def read_file(self, filename, opening) -> dict:
        """open fasta and return all the sequence in a dictionnary

        Args:
            filename (str): path to file

        Returns:
            dict: all the seq of the file
        """
        acid_amine = ['C', 'S', 'T', 'A', 'G', 'P', 'D', 'E', 'Q',
                      'N', 'H', 'R', 'K', 'M', 'I', 'L', 'V', 'W', 'Y', 'F']
        counter = 0
        with open(filename, f'{opening}') as self.file:
            for line in self.file:
                if not line.startswith(">") and line.strip().split() in acid_amine:
                    ValueError("Not a fasta file ! Please, verify your file")
                elif line.startswith(">"):
                    counter += 1
                    temp_header = line.strip().split('|')[0][1:]
                    self.dict_of_seq[temp_header] = ""
                else:
                    self.dict_of_seq[temp_header] += line.strip()
        self.alignement.dict_sequence_tree = self.dict_of_seq






if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Align and tree")
    parser.add_argument('file_name', metavar='file_name', type=str,
                    help='name of the fastafile')
    parser.add_argument('--output', metavar='output', type=str,
                    help='name of the output file', default='output.txt')
    args =  parser.parse_args()
    
    
    model = Model()
    
    # On lance les fonctions dans l'ordre nécessaire !!!
    model.read_file(args.file_name, "r")
    a = model.alignement.get_all_max_score(model.dict_of_seq)
    d = model.alignement.conserved_position()
    model.alignement.conserved_distance_matrix()
    model.alignement.nj_global()

    # Ecriture du fichier d'output    
    
    with open(args.output, "w") as file:
        
        # Mise en forme des résultats
        
        # Pour la matrice des scores :
        file.write("Matrice des scores : \n\n")
        for line in model.alignement.output_score_mat:
            file.write("[")
            for element in line:
                file.write(str(element))
                file.write("\t")
            file.write("]\n")
        file.write("\n\n\n")
        
        # Arbre guide newick :
        file.write("Arbre guide au format newick : \n\n")
        file.write(str(model.alignement.output_newick_upgma))
        file.write("\n\n\n")
        
        # Alignement multiple :
        file.write("\nAlignement multiple : \n\n")
        for i in range(len(model.alignement.output_clades_names_for_align)):
            file.write(str(model.alignement.output_clades_names_for_align[i])+
                       " : \t"+str(model.alignement.output_multiple_align[i])+"\n")
        file.write("\n\n")
        
        # Matrice positions conservées :
        file.write("\nMatrice des positions conservées : \n\n")
        for line in model.alignement.output_conserved_pos_mat:
            file.write("[")
            for element in line:
                file.write(str(element))
                file.write("\t")
            file.write("]\n")
        file.write("\n\n\n")
        
        # Final newick :
        file.write("Arbre newick arpès NJ: \n\n")
        file.write(str(model.alignement.output_newick_final))
        file.write("\n\n\n")
            
