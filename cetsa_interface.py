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


class MethodChooser(tk.Frame):
    
    NOT_SET = 0
    COMBO_T_TEST = 1
    NPARC = 2
    
    def __init__(self, master, variable):
        tk.Frame.__init__(self, master)
        
        self.variable = variable
        
        self.ita_radio = ttk.Radiobutton(self, 
                                         text="ITA", 
                                         variable=variable,
                                         value=self.COMBO_T_TEST)
        
        self.nparc_radio = ttk.Radiobutton(self,
                                           text="NPARC",
                                           variable=variable,
                                           value=self.NPARC)
        
        self.method_label = ttk.Label(self,
                                      text="Analysis method")
        
        self.method_label.grid(row=0, column=0)
        self.ita_radio.grid(row=1, column=0)
        self.nparc_radio.grid(row=2, column=0)
    
    def get(self):
        return self.variable.get()


window = tk.Tk(screenName="CETSA Pipeline")

approach = tk.IntVar(value=MethodChooser.NOT_SET)

method_choose = MethodChooser(window, approach)

candidate_filepath = tk.StringVar()

data_filepath = tk.StringVar()

candidate_ask = FilenameOpenChooser(window, 
                                    "Candidate", 
                                    candidate_filepath,
                                    cetsa_paths.get_candidates_filepath(False))

data_ask = FilenameOpenChooser(window,
                               "Data",
                               data_filepath,
                               cetsa_paths.get_data_filepath(False))

candidate_ask.pack()
data_ask.pack()
method_choose.pack()

window.mainloop()

print(candidate_filepath.get())
