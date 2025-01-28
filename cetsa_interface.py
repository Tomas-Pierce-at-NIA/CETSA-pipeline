# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 16:10:56 2024

@author: piercetf
"""

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox

import multiprocessing as mp
from pathlib import Path


import toml

import cetsa_paths
import CETSA_individual_temperature as ita
import cetsa2
import load_monocyte_cetsa_data as load




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
    
    COMBO_T_TEST = 1
    NPARC = 0
    
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


class ConditionList(tk.Frame):
    
    def __init__(self, master, titletxt, addtxt, condlist):
        super().__init__(master)
        
        self.idents = []
        self.next_idx = 0
        self.cond_list = condlist
        
        self.label = ttk.Label(self, text=titletxt)
        self.scrollbar = ttk.Scrollbar(self)
        self.listbox = ttk.Treeview(self, 
                                    yscrollcommand=self.scrollbar.set,
                                    show="tree")
        self.scrollbar.configure(command=self.listbox.yview)
        
        self.__newtext = tk.StringVar()
        
        self.addlabel = ttk.Label(self, text=addtxt)
        self.textbox = ttk.Entry(self, textvariable=self.__newtext)
        
        
        self.add_button = ttk.Button(self, 
                                     text='Add condition',
                                     command=self.add_condition)
        
        self.remove_button = ttk.Button(self,
                                        text="Remove condition",
                                        command=self.remove_condition)
        
        self.label.grid(row=0, column=0, columnspan=4, sticky='W')
        self.scrollbar.grid(row=1, column=4, rowspan=4)
        self.listbox.grid(row=1, column=0, columnspan=4, rowspan=4)
        self.addlabel.grid(row=5, column=0, columnspan=4, sticky='W')
        self.textbox.grid(row=6, column=0, columnspan=4, sticky='W')
        self.add_button.grid(row=7, column=0, columnspan=2)
        self.remove_button.grid(row=7, column=2, columnspan=2)
        
        
        self.__condition_add_initial()
    
    def __condition_add_initial(self):
        for item in self.cond_list:
            identifier = self.listbox.insert("", "end", self.next_idx, text=item)
            self.idents.append(identifier)
            self.listbox.see(identifier)
            self.next_idx += 1
    
    
    def add_condition(self):
        text = self.__newtext.get()
        identifier = self.listbox.insert("", "end", self.next_idx, text=text)
        self.idents.append(identifier)
        self.listbox.see(identifier)
        self.next_idx += 1
        self.cond_list.append(text)
    
    
    def remove_condition(self):
        if len(self.idents) > 0:
            last = self.idents.pop()
            self.listbox.delete(last)
            self.cond_list.pop()
        else:
            messagebox.showerror(title="Error",
                                 message="No condition to remove")


class ConditionEntry(tk.Frame):
    
    def __init__(self, master, variable, titletxt):
        super().__init__(master)
        
        self.variable = variable
        
        self.label = ttk.Label(self, text=titletxt)
        self.entry = ttk.Entry(self, textvariable=self.variable)
        
        self.label.grid(row=0, column=0, sticky='W')
        self.entry.grid(row=1, column=0, sticky='W')
    
    def get(self):
        return self.variable.get()


class Config:
    
    def __init__(self):
        self.__prev_vcontrol = None
        self.__prev_nscontrols = None
        self.__prev_outdir = None
        self.__prev_logdir = None
        self.__prev_cachedata = None
        self.__prev_cachecandidate = None
        self.__prev_unc_data = None
        self.__prev_unc_candidate = None
        self.__use_cached = None
        self.load_prev()
        
        self.vcontrol = tk.StringVar()
        self.ns_controls = []
        self.outdir = tk.StringVar()
        self.logdir = tk.StringVar()
        self.data = tk.StringVar()
        self.candidate = tk.StringVar()
        self.usecache = tk.BooleanVar()
        self.use_prev()
        #self.use_cached()
        
        self.method = tk.IntVar()
    
    
    def display(self):
        print("vehicle", self.vcontrol.get())
        print("nonsenolytics", self.ns_controls)
        print("outdir", self.outdir.get())
        print("logdir", self.logdir.get())
        print("data", self.data.get())
        print("candidate", self.candidate.get())
        print("use cache", self.usecache.get())
        print("method", self.method.get())
    
    
    @property
    def default_logdir(self):
        return self.__prev_logdir
    
    
    @property
    def default_outdir(self):
        return self.__prev_outdir
    
    @property
    def default_data(self):
        return self.__prev_cachedata
    
    
    @property
    def default_candidate(self):
        return self.__prev_cachecandidate
    
    
    def ready_config(self):
        path = cetsa_paths.paramfilename()
        use_cache = self.usecache.get()
        data = {
            'title' : "CETSA pipeline configuration file",
            'controls': {
                'vehicle': self.vcontrol.get(),
                'nonsenolytic': self.ns_controls
                },
            'outputs': {
                'outdir': self.outdir.get(),
                'logdir': self.logdir.get()
                },
            'inputs': {'use_cache': use_cache}
            }
        if use_cache:
            data['inputs']['cached'] = {
                'data': self.data.get(),
                'candidates': self.candidate.get()
                }
            data['inputs']['uncached'] = {
                'data': self.data.get(),
                'candidates': self.candidate.get()
                }
        
        else:
            data['inputs']['cached'] = {
                'data': self.data.get(),
                'candidates': self.candidate.get()
                }
            data['inputs']['uncached'] = {
                'data': self.data.get(),
                'candidates': self.candidate.get()
                }
        
        with open(path, 'w') as confighandle:
            toml.dump(data, confighandle)
    
    
    
    def use_prev(self):
        self.vcontrol.set(self.__prev_vcontrol)
        self.ns_controls = self.__prev_nscontrols[:]
        self.outdir.set(self.__prev_outdir)
        self.logdir.set(self.__prev_logdir)
        self.data.set(self.__prev_unc_data)
        self.candidate.set(self.__prev_unc_candidate)
        self.usecache.set(self.__use_cached)
    
    
    def use_cached(self):
        self.data.set(self.__prev_cachedata)
        self.candidate.set(self.__prev_cachecandidate)
        self.usecache.set(True)
    
    
    def use_uncached(self):
        self.data.set(self.__prev_unc_data)
        self.candidate.set(self.__prev_cachecandidate)
        self.usecache.set(False)
    
    
    def load_prev(self):
        oldparams = cetsa_paths.loadparams()
        self.__prev_vcontrol = oldparams['controls']['vehicle']
        self.__prev_nscontrols = oldparams['controls']['nonsenolytic']
        self.__prev_outdir = oldparams['outputs']['outdir']
        self.__prev_logdir = oldparams['outputs']['logdir']
        self.__prev_cachedata = oldparams['inputs']['cached']['data']
        self.__prev_cachecandidate = oldparams['inputs']['cached']['candidates']
        self.__prev_unc_data = oldparams['inputs']['uncached']['data']
        self.__prev_unc_candidate = oldparams['inputs']['uncached']['candidates']
        self.__use_cached = oldparams['inputs']['use_cache']
        


def ita_wrapper():
    candidatepath = cetsa_paths.candidate_filename()
    datapath = cetsa_paths.data_filename()
    outdir = cetsa_paths.get_outdir()
    data, can = load.prepare_data(display=False, 
                      data_path=datapath, 
                      candidate_path=candidatepath)
    ita.run_analysis(data, can, outdir)


    
def ui():
    window = tk.Tk()
    
    configdata = Config()
    
    
    cachecheck = ttk.Checkbutton(window, variable=configdata.usecache, text='cached')
    datachoice = FilenameOpenChooser(window, 
                                     "data", 
                                     configdata.data,
                                     defaultpath=configdata.default_data)
    candidatechoice = FilenameOpenChooser(window,
                                          "candidate",
                                          configdata.candidate,
                                          defaultpath=configdata.default_candidate)
    
    outdirchoice = FolderOpenChooser(window,
                                     "output",
                                     configdata.outdir,
                                     defaultpath=configdata.default_outdir)
    
    logdirchoice = FolderOpenChooser(window,
                                     "log",
                                     configdata.logdir,
                                     defaultpath=configdata.default_logdir)
    
    vehicleentry = ConditionEntry(window, 
                                  configdata.vcontrol, 
                                  "Vehicle control")
    
    nonsenolyticlist = ConditionList(window,
                                     'Nonsenolytic controls',
                                     'Nonsenolytic',
                                     configdata.ns_controls)
    
    titlelabel = ttk.Label(window,
                           text="SenoCETSA graphical user interface")
    
    methodchoice = MethodChooser(window, configdata.method)
    
    displaybutton = ttk.Button(window,
                               text="Display parameters",
                               command=configdata.display)
    
    quitbutton = ttk.Button(window,
                            text="Quit",
                            command=window.destroy)
    
    runbutton = ttk.Button(window,
                           text="Run",
                           command=messagebox.showinfo(
                               title="Not implemented",
                               message="Running has not been implemented")
                           )
    
    
    #cachecheck.grid(row=0, column=0)
    titlelabel.grid(row=0, column=0, columnspan=5)
    quitbutton.grid(row=1, column=1)
    cachecheck.grid(row=1, column=2)
    datachoice.grid(row=2,column=1, columnspan=2)
    candidatechoice.grid(row=3, column=1, columnspan=2)
    outdirchoice.grid(row=4, column=1, columnspan=2)
    logdirchoice.grid(row=5, column=1, columnspan=2)
    
    methodchoice.grid(row=6, column=1)
    
    vehicleentry.grid(row=1, column=3, columnspan=2, sticky='W')
    nonsenolyticlist.grid(row=2, column=3, columnspan=2, rowspan=6)
    
    displaybutton.grid(row=7, column=1)
    runbutton.grid(row=7, column=2)
    
    window.title("SenoCETSA")
    window.mainloop()


if __name__ == '__main__':
    ui()

