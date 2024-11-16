# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 16:10:56 2024

@author: piercetf
"""

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog

import cetsa_paths

class FilenameOpenChooser(tk.Frame):
    
    def __init__(self, master, name: str, textvariable, defaultpath=None):
        tk.Frame.__init__(self, master)
        
        self.textvariable = textvariable
        self.defaultpath = defaultpath
        self.label = ttk.Label(self,
                               text=f"Input {name} filename")
        self.entry = ttk.Entry(self, 
                               textvariable=textvariable)
        self.button = ttk.Button(self, 
                                 command=self.choose, 
                                 text=f"Choose {name} filename")
        self.reset = ttk.Button(self,
                                 command=self.default,
                                 text=f"Default {name} filename")
        self.label.grid(row=0, column=0)
        self.entry.grid(row=0, column=1)
        self.button.grid(row=1, column=0)
        self.reset.grid(row=1, column=1)
        
    
    def choose(self):
        filepath = filedialog.askopenfilename()
        self.textvariable.set(filepath)
        
    def default(self):
        if self.defaultpath is not None:
            self.textvariable.set(self.defaultpath)
        else:
            self.textvariable.set("")
    
    def get(self):
        return self.textvariable.get()

NOT_SET = 0
COMBO_T_TEST = 1
NPARC = 2


window = tk.Tk(screenName="CETSA Pipeline")

approach = tk.IntVar(value=NOT_SET)

candidate_filepath = tk.StringVar()


methodlab = ttk.Label(window, text="Analysis method")
itabut = ttk.Radiobutton(master=window, variable=approach, text="ITA", value=COMBO_T_TEST)
nparcbut = ttk.Radiobutton(master=window, variable=approach, text="NPARC", value=NPARC)

candidate_ask = FilenameOpenChooser(window, 
                                    "Candidate", 
                                    candidate_filepath,
                                    cetsa_paths.get_candidates_filepath(False))

candidate_ask.pack()
methodlab.pack()
itabut.pack()
nparcbut.pack()

window.mainloop()

print(candidate_filepath.get())
