# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 16:10:56 2024

@author: piercetf
"""

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog

import cetsa_paths


class FsChooser(tk.Frame):
    
    def __init__(self, master, name: str, textvariable, command, kind="filename", defaultpath=None):
        tk.Frame.__init__(self, master)
        self.command = command
        self.textvariable = textvariable
        self.defaultpath = defaultpath
        self.label = ttk.Label(self,
                               text=f"Input {name} {kind}")
        self.entry = ttk.Entry(self, 
                               textvariable=textvariable)
        self.button = ttk.Button(self, 
                                 command=self.choose, 
                                 text=f"Choose {name} {kind}")
        self.reset = ttk.Button(self,
                                 command=self.default,
                                 text=f"Default {name} {kind}")
        self.label.grid(row=0, column=0)
        self.entry.grid(row=0, column=1)
        self.button.grid(row=1, column=0)
        self.reset.grid(row=1, column=1)
        
    
    def choose(self):
        filepath = self.command()
        self.textvariable.set(filepath)
        
    def default(self):
        if self.defaultpath is not None:
            self.textvariable.set(self.defaultpath)
        else:
            self.textvariable.set("")
    
    def get(self):
        return self.textvariable.get()


class FolderOpenChooser(FsChooser):
    
    def __init__(self, 
                 master, 
                 name: str, 
                 textvariable, 
                 command=filedialog.askdirectory, 
                 kind="folder", 
                 defaultpath=None):
        super().__init__(master, name, textvariable, command, kind, defaultpath)





class FilenameOpenChooser(FsChooser):
    
    def __init__(self, 
                 master, 
                 name: str, 
                 textvariable, 
                 command=filedialog.askopenfilename, 
                 kind="file", 
                 defaultpath=None):
        super().__init__(master, name, textvariable, command, kind, defaultpath)
        


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

window.title("CETSA Pipelines")

approach = tk.IntVar(value=MethodChooser.NOT_SET)

method_choose = MethodChooser(window, approach)

candidate_filepath = tk.StringVar()

data_filepath = tk.StringVar()

outdir_filepath = tk.StringVar()

candidate_ask = FilenameOpenChooser(window, 
                                    "Candidate", 
                                    candidate_filepath,
                                    defaultpath=cetsa_paths.get_candidates_filepath(False))

data_ask = FilenameOpenChooser(window,
                               "Data",
                               data_filepath,
                               defaultpath=cetsa_paths.get_data_filepath(False))

out_ask = FolderOpenChooser(window,
                            "Output",
                            outdir_filepath,
                            defaultpath=str(cetsa_paths.get_outdir())
                            )

quit_but = ttk.Button(window, command=window.destroy, text="Quit")

candidate_ask.pack()
data_ask.pack()
out_ask.pack()
method_choose.pack()
quit_but.pack()

window.mainloop()

print(candidate_filepath.get())
print(data_filepath.get())
print(approach.get())
