from arbre import *
from constants import *
from alignement import *


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
    model = Model()
    # On lance les fonctions dans l'ordre nécessaire !!!
    model.read_file("opsines_juste4.fasta.txt", "r")
    a = model.alignement.get_all_max_score(model.dict_of_seq)
    b = model.alignement.upgma.tree_with_upgma(a)
    print(b)
    model.alignement.multiple_alignement()
