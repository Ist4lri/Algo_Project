from model import *
from vue import Application


class Controller():
    def __init__(self, filename="opsines_juste4.fasta.txt"):
        """Init of the Controller"""
        self.model = Model()
        self.model.read_file(filename, "r")
        self.view = Application(self)
        self.view.view_window()

    def open_other_file(self, filename):
        self.model.file_close()
        self.model.read_file(filename, opening="r")


c = Controller()
