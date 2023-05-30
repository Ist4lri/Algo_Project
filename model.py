from arbre import *
from constants import *
from alignement import *
import argparse
from Save_to_output import writing_output_file


def is_fasta(fastafile):  # takes in argument the name of the file
    """This function control if the file starts with a '>' used in fasta files headers."""
    if fastafile.split(".")[-1] in ["fasta", "faa", "fa"]:  # extension is OK
        with open(fastafile, "r", encoding='UTF-8') as file:
            first_line = file.read(1)  # take the first line of the file
            if first_line.split()[0] == ">":  # first char of first line
                return True  # this looks like a fasta file
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

    def app_guidelines(self):
        """Launching all the Script, in right order"""
        self.alignement.get_all_max_score(self.dict_of_seq)
        self.alignement.conserved_position()
        self.alignement.conserved_distance_matrix()
        self.alignement.nj_global()

    def show_matrixes(self):
        """Return all Matrix"""
        return self.alignement.output_score_mat, self.alignement.output_conserved_pos_mat

    def show_trees(self):
        """Return all Trees"""
        return self.alignement.output_newick_upgma, self.alignement.output_newick_final

    def show_alignements(self):
        """Return all Alignements"""
        return self.alignement.output_clades_names_for_align, self.alignement.output_multiple_align


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Here's a program for making multiple alignments, calculating scores, and building trees based on these scores.")
    parser.add_argument('-f', '--file_name', metavar='Path of Fasta File', type=str,
                        help='name of the fastafile')
    parser.add_argument('-o', '--output', metavar='Output file path', type=str,
                        help='name of the output file', default='output.txt', required=False)
    args = parser.parse_args()

    model = Model()

    # On lance les fonctions dans l'ordre nécessaire !!!
    model.read_file(args.file_name, "r")
    model.app_guidelines()

    # Ecriture du fichier d'output
    writing_output_file(args.output, model)
