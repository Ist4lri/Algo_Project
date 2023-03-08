from arbre import *
from constants import *
from alignement import *


class Model():
    """Model Class"""

    def __init__(self):
        """Init of the class"""
        self.dict_of_seq = {}
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
        with open(filename, f'{opening}') as self.file:
            for line in self.file:
                if not line.startswith(">") and line.strip().split() in acid_amine:
                    ValueError("Not a fasta file ! Please, verify your file")
                elif line.startswith(">"):
                    temp_header = line.strip()
                    self.dict_of_seq[line.strip()] = ""
                else:
                    self.dict_of_seq[temp_header] += line.strip()
        print(self.alignement.get_all_max_score(self.dict_of_seq))

    def send_dict(self) :
        self.alignement.get_all_max_score(self.dict_of_seq)
        
if __name__ == "__main__" :
    model = Model()
    model.read_file("opsines_juste4.fasta.txt", "r")
    model.alignement.get_all_max_score(model.dict_of_seq)
    print(model.alignement.upgma.tree_with_upgma(model.alignement.upgma, model.alignement.matrice_distance, model.alignement.clades_names))
    