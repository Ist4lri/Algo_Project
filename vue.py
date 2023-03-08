"""For the Visual part of the app"""
from tkinter import Listbox, Button, Tk, Label, Entry, LabelFrame
from tkinter.messagebox import showinfo, showerror
from tkinter import filedialog


class Application(Tk):
    """The Application himself, the design of the App"""

    def __init__(self, controller) -> None:
        """Initiation of class

        Args:
            controller (object): The controller of the app
        """
        Tk.__init__(self)
        self.initiate_app()
        self.controller = controller

    def initiate_app(self):
        """Initiation of the graphic phase of app"""
        self.subLabel = LabelFrame(self)
        self.explanation = Label(
            self.subLabel, text="Voici un logiciel pour pouvoir faire des alignements multiples, calculer les score, pouvoir constuire des arbres en fonction de ces score. \nDifférentes méthodes ont été utilisées: (liste des méthodes.) \nVeuillez charger un Jeu de donnée avant de  lancer toutes actions.\nPour information, un jeu de donnée par défaut est déjà enregistré.").pack()
        self.charge_file_b = Button(
            self.subLabel, text="Add a file", command=self.open_file).pack()
        self.alignement_b = Button(
            self.subLabel, text="Show Alignement", command=self.show_alignement).pack()
        self.tree_b = Button(
            self.subLabel, text="Show Tree", command=self.show_tree).pack()
        self.quit_b = Button(self.subLabel, text="Quit",
                             command=self.destroy).pack()

        self.subLabel.pack(fill="both", expand="yes", padx=30, pady=20)

    def open_file(self):
        """To Open another file than the default fasta file"""
        filename = filedialog.askopenfilename(
            title="Chosse File", filetypes=(("fasta file", "faa"), ("fasta file", "fasta"), ("text file", "txt")))
        self.controller.open_other_file(filename)

    def show_alignement(self):
        print("bite")

    def show_tree(self):
        print("chatte")

    def view_window(self):
        """Launch the Window"""
        self.title("Alignement/Arbre app")
        self.mainloop()
