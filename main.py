from model import *
from vue import Application


class Controller():
    def __init__(self, filename="opsines.fasta.txt"):
        """Init of the Controller"""
        self.model = Model()
        self.model.read_file(filename, "r")
        self.view = Application(self)
        self.view.view_window()

    def open_other_file(self, filename):
        self.model.read_file(filename, opening="r")

    def launch_app(self):
        self.model.app_guidelines()

    def show_all_matrix(self):
        mat_1, mat_2 = self.model.show_matrixes()
        self.view.display_result(
            "Matrice des scores", mat_1, "Matrice des positions conservées", mat_2)

    def show_all_trees(self):
        tree_1, tree_2 = self.model.show_trees()
        self.view.display_result(
            "Arbre guide",
            tree_1,
            "Arbre final après NJ",
            tree_2)

    def show_all_alignements(self):
        MSA, MSA_content = self.model.show_alignements()
        output = []
        for i in range(len(MSA)):
            output.append(str(MSA[i]) +
                          " : \t"+str(MSA_content[i])+"\n")
        self.view.display_result("Alignement Multiple", "".join(output))


if __name__ == "__main__":
    c = Controller()
