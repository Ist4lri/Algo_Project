from model import *
from vue import Application


class Controller():
    def __init__(self, filename="opsines_juste4.fasta.txt"):
        """Init of the Controller"""
        self.model = Model()
        print(self.model.read_file(filename, "r"))
        self.view = Application(self)
        self.view.view_window()


c = Controller()
