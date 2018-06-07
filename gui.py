import tkinter as tk
from tkinter import ttk
import subprocess
import compare

def run_script():
    subprocess.run(["python", "run_compare.py"])

def run_compare():
    compare.start("compare_output1.xml", "leptolyngbya.fasta")

root = tk.Tk()
root.title("SeqComapre")

main_frame = ttk.Frame(root, padding="3 3 12 12")

main_frame.grid(column=0, row=0, sticky=(tk.N, tk.W, tk.E, tk.S))
main_frame.columnconfigure(0, weight=1)
main_frame.rowconfigure(0, weight=1)

button = ttk.Button(main_frame, text='Go', command=run_script)
#button.grid(column=3, row=3, sticky=tk.W)

button_compare = ttk.Button(main_frame, text='compare', command=run_compare)
#button_compare.grid(column=3, row=3, sticky=tk.W)

for child in main_frame.winfo_children(): child.grid_configure(padx=5, pady=5)

root.mainloop()