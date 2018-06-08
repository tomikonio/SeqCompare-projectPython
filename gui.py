import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import subprocess
import compare
import os


class Gui(ttk.Frame):
    def __init__(self, root, ):
        super().__init__(master=root, padding="3 3 12 12")

        self.root = root
        self.file_dict = {}
        self.file_name_labels = []

        self.pack()
        self.create_widgets()
        for child in self.winfo_children(): child.grid_configure(padx=5, pady=5)

    def run_script(self):
        subprocess.run(["python", "run_compare.py"])

    def run_compare(self):
        compare.start("compare_output1.xml", "leptolyngbya.fasta")

    def select_folder(self):
        folder_path = filedialog.askdirectory()
        print(folder_path)

        # Todo check if a folder was selected

        self.file_dict.clear()
        for file in os.scandir(folder_path):
            if file.name.endswith('.fasta'):
                self.file_dict[file.name] = file.path
        print(self.file_dict)
        self.create_file_labels()

    def create_file_labels(self):
        for label in self.file_name_labels:
            label.destroy()
        self.file_name_labels.clear()

        for file_name in self.file_dict:
            label = ttk.Label(self, text=file_name)
            label.grid()
            self.file_name_labels.append(label)

    def create_widgets(self):
        # button = ttk.Button(self, text='Go', command=self.run_script)
        # button.grid(column=3, row=3, sticky=tk.W)

        # button_compare = ttk.Button(self, text='compare', command=self.run_compare)
        # button_compare.grid(column=3, row=3, sticky=tk.W)

        label = ttk.Label(self, text='Select a folder:').grid(column=1,row=2)
        button_folder = ttk.Button(self, text='Browse...', command=self.select_folder).grid(column=2, row=2)


root = tk.Tk()
root.title("SeqCompare")
app = Gui(root)
app.mainloop()

# file_dict = {}
#
# root = tk.Tk()
# root.title("SeqComapre")
#
# main_frame = ttk.Frame(root, padding="3 3 12 12")
#
# main_frame.grid(column=0, row=0, sticky=(tk.N, tk.W, tk.E, tk.S))
# main_frame.columnconfigure(0, weight=1)
# main_frame.rowconfigure(0, weight=1)
#
# button = ttk.Button(main_frame, text='Go', command=run_script)
# # button.grid(column=3, row=3, sticky=tk.W)
#
# button_compare = ttk.Button(main_frame, text='compare', command=run_compare)
# # button_compare.grid(column=3, row=3, sticky=tk.W)
#
#
# label = ttk.Label(main_frame, text='Select a folder:')
# button_folder = ttk.Button(main_frame, text='Browse...', command=select_folder)
#
# for child in main_frame.winfo_children(): child.grid_configure(padx=5, pady=5)
#
# root.mainloop()
