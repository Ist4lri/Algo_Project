"""For the Visual part of the app"""
from tkinter import Button, Tk, Label, LabelFrame, Toplevel, StringVar, Entry
from tkinter.messagebox import showinfo
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
        self.win_width = 0

    def initiate_app(self):
        """Initiation of the graphic phase of app"""

        self.win_width = self.winfo_screenwidth()

        self.subLabel = LabelFrame(self)
        Label(
            self.subLabel, text="Here's a program for making multiple alignments, calculating scores, and building trees based on these scores. \nVarious methods were used: (list of methods.) \nPlease load a dataset before launching any action. \nFor information, a default dataset is already saved.").pack()
        self.charge_file_b = Button(
            self.subLabel, text="Add a file", command=self.open_file).pack()
        Label(self.subLabel, text="\n").pack()
        self.launch_b = Button(
            self.subLabel, text="Launch the MSA", command=self.launch).pack()
        Label(self.subLabel, text="\n").pack()
        self.quit_b = Button(self.subLabel, text="Quit",
                             command=self.destroy).pack()

        self.subLabel.pack(fill="both", expand="yes", padx=30, pady=20)

    def open_file(self):
        """To Open another file than the default fasta file"""
        filename = filedialog.askopenfilename(
            title="Chosse File", filetypes=(("fasta file", ".faa"), ("fasta file", ".fasta"), ("text file", ".txt")))
        self.controller.open_other_file(filename)

    def launch(self):
        """For launch the MSA script"""
        showinfo("Information",
                 "MSA running take some minutes.\nPlease, wait for the other Window.")
        self.controller.launch_app()
        window2 = Toplevel()
        Button(
            window2, text="Show Matrices", command=self.show_all_matrix).pack()
        Label(window2, text="\n").pack()
        Button(
            window2, text="Show Alignements", command=self.show_all_alignement).pack()
        Label(window2, text="\n").pack()
        Button(
            window2, text="Show Trees", command=self.show_all_tree).pack()
        Label(window2, text="\n").pack()
        Button(window2, text="Quit", command=window2.destroy).pack()

    def show_all_alignement(self):
        """For showing the MSA"""
        self.controller.show_all_alignements()

    def show_all_tree(self):
        """For showing the Trees"""
        self.controller.show_all_trees()

    def show_all_matrix(self):
        """For showing all Matrix"""
        self.controller.show_all_matrix()

    def display_result(self, first_name="", first_to_display="", second_name="", second_to_display=""):
        """To display The result on a another Window

        Args:
            first_name (str, optional): Name of first result. Defaults to "".
            first_to_display (str, optional): The content of first result. Defaults to "".
            second_name (str, optional): Name of second result. Defaults to "".
            second_to_display (str, optional): The content of second result. Defaults to "".
        """
        value1 = StringVar()
        value2 = StringVar()
        value1.set(first_to_display)
        value2.set(second_to_display)
        another_window = Toplevel()
        Label(another_window, text=first_name).pack()
        Entry(another_window, textvariable=value1,
              cursor="target", width=100).pack()
        Label(another_window, text="\n\n").pack()
        if second_to_display != "":
            Label(another_window, text=second_name).pack()
            Entry(another_window, textvariable=value2,
                  cursor="target", width=100).pack()
            Label(another_window, text="\n\n").pack()
        Button(another_window, text="Quit this Window",
               command=another_window.destroy).pack()

    def view_window(self):
        """Launch the Window"""
        self.title("Alignement/Arbre app")
        self.mainloop()
