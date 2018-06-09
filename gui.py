import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
import subprocess
import compare
import run_compare
import os
from collections import OrderedDict


# TODO make an option to select in which order the files will be compared
# TODO check if the numbering is working

class Gui(ttk.Frame):
    def __init__(self, root, ):
        super().__init__(master=root, padding="3 3 12 12")

        self.root = root
        self.file_dict = {}  # store file-name: file-path pairs
        self.file_name_labels = []  # store label objects that are associated with secondary file names
        self.file_combos = {}   # store combobox objects for secondary files
        self.folder_path = ""
        self.primary_combo = ""
        self.primary_file = ""
        # self.label_frame = ""
        self.number_of_files = 0

        self.pack()
        self.create_widgets()
        self.new_frame = ""
        self.padd()

    def padd(self):
        for child in self.winfo_children(): child.grid_configure(padx=5, pady=5)

    def check_combo_input(self):
        for file_name in self.file_combos:
            # print(self.file_combos[file_name].get())
            if self.file_combos[file_name][1].current() == -1:
                messagebox.showinfo("Error", "Please select match/not match for all files")
                return False
            elif self.file_combos[file_name][0].current() == -1:
                messagebox.showinfo("Error", "Please select a number for all files")
                return False

        return True

    def run_script(self):
        # os.chdir(self.folder_path)
        if self.check_combo_input():
            secondary_files = OrderedDict()
            for i in range(1, self.number_of_files):
                for file_name in self.file_combos:
                    #print(self.file_combos[file_name].get())
                    # if self.file_combos[file_name][1].current() == -1:
                    #     messagebox.showinfo("Error", "Please select match/not match for all files")
                    # elif self.file_combos[file_name][0].current() == -1:
                    #     messagebox.showinfo("Error", "Please select a number for all files")
                    if self.file_combos[file_name][0].get() == str(i):
                        secondary_files[file_name] = "m" if self.file_combos[file_name][1].get() == "match" else "nm"
            print(secondary_files)
            run_compare.start(self.primary_file,secondary_files, self.folder_path)
        # subprocess.run(["python", "run_compare.py"])

    # def run_compare(self):
    #     compare.start("compare_output1.xml", "leptolyngbya.fasta")

    def select_folder(self):
        new_folder_path = filedialog.askdirectory()
        print(new_folder_path)

        if new_folder_path:
            self.folder_path = new_folder_path
            self.file_dict.clear()
            for file in os.scandir(new_folder_path):
                if file.name.endswith('.fasta'):
                    self.file_dict[file.name] = file.path
            if len(self.file_dict.keys()) < 2:
                messagebox.showinfo("Error", "Please select a folder with at least two .fasta files")
            else:
                print(self.file_dict)
                self.number_of_files = len(self.file_dict.keys())
                ttk.Label(self, text="The current folder is {}".format(self.folder_path)).grid(column=3, row=2)
                self.choose_primary()
                # self.create_file_labels()

    def primary_selected(self, event):
        self.primary_file = self.primary_combo.get()
        print(self.primary_file)
        #ttk.Separator(self, orient=tk.HORIZONTAL)
        self.create_file_labels()


    def choose_primary(self):
        self.destroy_things()
        list = []
        for file_name in self.file_dict:
            list.append(file_name)
        ttk.Label(self, text="Select a Primary file").grid(column=1, row=5)
        self.primary_combo = ttk.Combobox(self, values=list, state='readonly')
        self.primary_combo.bind('<<ComboboxSelected>>', self.primary_selected)
        self.primary_combo.grid(column=2, row=5)
        self.padd()

    def combo_number_change(self, event):
        print(event.widget)
        #
        # chosen_number = event.widget.get()
        # num_list = []
        #
        # for file_name in self.file_combos:
        #     number_combo = self.file_combos[file_name][0]
        #     if number_combo != event.widget:
        #         number_combo['values'] =

    def destroy_things(self):
        if self.new_frame != "":
            self.new_frame.destroy()
        for file_name in self.file_combos:
            self.file_combos[file_name][0].destroy()
            self.file_combos[file_name][1].destroy()
        self.file_combos.clear()

        for label in self.file_name_labels:
            label.destroy()
        self.file_name_labels.clear()

    def create_file_labels(self):
        row = 7
        column = 1

        self.destroy_things()
        # if self.new_frame != "":
        #     self.new_frame.destroy()
        # for file_name in self.file_combos:
        #     self.file_combos[file_name][0].destroy()
        #     self.file_combos[file_name][1].destroy()
        # self.file_combos.clear()
        #
        # for label in self.file_name_labels:
        #     label.destroy()
        # self.file_name_labels.clear()

        self.new_frame = ttk.Labelframe(self, text="Fasta files")
        self.new_frame.grid(column=1,row=6, columnspan=3)
        # self.label_frame = ttk.Labelframe(self, text="FASTA files")
        # self.label_frame.grid(column=1, row=6)

        for file_name in self.file_dict:
            label = ttk.Label(self.new_frame, text=file_name)
            label.grid(column=column, row=row)
            if file_name == self.primary_file:
                ttk.Label(self.new_frame, text="Primary").grid(column=column + 1, row=row)
            else:
                combo = ttk.Combobox(self.new_frame, values=["match", "not match"], state='readonly')
                combo.grid(column=column + 1, row=row)

                number_list = list(range(1, self.number_of_files))
                number_combo = ttk.Combobox(self.new_frame, values=number_list, state='readonly')
                number_combo.bind('<<ComboboxSelected>>', self.combo_number_change)
                number_combo.grid(column=column + 2, row=row)

                self.file_combos[file_name] = [number_combo, combo]

            row += 1

            self.file_name_labels.append(label)
        button = ttk.Button(self.new_frame, text='Go', command=self.run_script)
        button.grid(column=column, row=row)
        for child in self.new_frame.winfo_children(): child.grid_configure(padx=5, pady=5)
        self.padd()

    def create_widgets(self):
        # button = ttk.Button(self, text='Go', command=self.run_script)
        # button.grid(column=3, row=3, sticky=tk.W)

        # button_compare = ttk.Button(self, text='compare', command=self.run_compare)
        # button_compare.grid(column=3, row=3, sticky=tk.W)

        label = ttk.Label(self, text='Select a folder:').grid(column=1, row=2)
        button_folder = ttk.Button(self, text='Browse...', command=self.select_folder).grid(column=2, row=2)
        ttk.Separator(self, orient=tk.HORIZONTAL)


root = tk.Tk()
root.title("SeqCompare")
root.minsize(800, 480)
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
