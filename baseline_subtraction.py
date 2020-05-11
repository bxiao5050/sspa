from tkinter import *
from tkinter import ttk
from tkinter.filedialog import askopenfilenames, askopenfilename, askdirectory, asksaveasfilename
import pandas as pd
import matplotlib
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure
import os
import numpy as np

from scipy import sparse
from scipy.sparse.linalg import spsolve
from tkinter import messagebox
import peakutils
import subprocess
from datetime import *

from rangeDragg import RangeDrag

class Baseline_subtraction(Frame):
    def __init__(self, master):
        super().__init__(master)
        self.master = master

        self.baseline_choose = '' #check which baseline method is used
                                    # 'auto', 'ALS'
        self.line = None
        self.line_baseline = None
        self.line_corrected = None

        self.alsPara = None # ALS parameter
        self.lam_l = None
        self.p_l = None

        self.expData = pd.DataFrame()
        leftFrame = Frame(self) # contain combobox+canvas
        leftFrame.pack(side = 'left', fill = 'both', expand = True)
        lf1 = LabelFrame(leftFrame, text = 'select data')
        self.cb= ttk.Combobox(lf1, width = 60)
        self.cb.pack()

        f = Figure(figsize=(5,5), dpi=100)
        self.ax = f.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(f, master = leftFrame)
        self.ax.set_xlabel('2Theta')
        self.ax.set_ylabel('intensity')

        toolbar = NavigationToolbar2Tk(self.canvas, leftFrame)
        toolbar.update()



        self.menubar()

    def menubar(self):
        self.menubar = Menu(self.master)

        #import
        importmanu = Menu(self.menubar)
        self.menubar.add_cascade(label="Import", menu=importmanu)
        importmanu.add_command(label="Import multiple .xy files", command=self.import_xy_files)
        importmanu.add_command(label="Import single .csv file", command=self.import_csv)

        #baseline
        baselinemanu = Menu(self.menubar)
        self.menubar.add_cascade(label="Baseline", menu=baselinemanu)
        baselinemanu.add_command(label="baseline auto", command=lambda method = 'auto': self.baseline(method))
        baselinemanu.add_command(label="Asymmetric Least Squares Smoothing", command=lambda method = 'ALS': self.baseline(method))

        self.menubar.entryconfig("Baseline", state = 'disabled')

        #export
        exportmanu = Menu(self.menubar)
        self.menubar.add_cascade(label="Export", menu=exportmanu)
        exportmanu.add_command(label="Export multiple .xy files", command=self.export_xy_files)
        exportmanu.add_command(label="Export single .csv file", command=self.export_csv)
        self.menubar.entryconfig("Export", state = 'disabled')

        #document
        documentmanu = Menu(self.menubar)
        self.menubar.add_cascade(label="Document", menu=documentmanu)
        documentmanu.add_command(label="Open references", command=self.openReference)


        # display the menu
        self.master.config(menu=self.menubar)

    def openReference(self):
        subprocess.Popen('443199618.pdf', shell = True)

    def on_updateXrange(self):
        self.drag.set_Xrange(left = float(self.xLeft_e.get()), right = float(self.xRight_e.get()))
        self.baseline(self.baseline_choose)

    def import_csv(self):
        path = askopenfilename(parent=self ,title='Choose a csv file', filetypes = (("csv files","*.csv"),("all files","*.*")))
        #
        if len(path) != 0:
            self.title = path
            basename = os.path.basename(path)
            filebase = os.path.splitext(basename)[1]
            # self.config(text = os.path.splitext(basename)[0])
        if len(path) != 0 and filebase == '.csv':
            self.expData = pd.read_csv(path, header = 0)

        self.polular_combobox(self.expData)

    def import_xy_files(self):
        files = askopenfilenames(parent=self,title='Choose .xy files', filetypes = (("xy files","*.xy"),("all files","*.*")))
        sortedfile = sorted([f for f in self.tk.splitlist(files)])

        for index, f in enumerate(sortedfile):
            df = pd.read_csv(f, header = 0, sep = ' ')
            if index == 0:
                self.expData.insert(0, 'Angle', df.iloc[:,0])
                self.expData.insert(1, os.path.splitext(os.path.basename(f))[0], df.iloc[:,1])
            else:
                self.expData.insert(index + 1, os.path.splitext(os.path.basename(f))[0], df.iloc[:,1])
        self.polular_combobox(self.expData)



    def baseline(self, method):
        if self.line_baseline is not None and len(self.line_baseline)>0:
            self.line_baseline.pop(0).remove()
            self.line_corrected.pop(0).remove()

        self.baseline_choose = method
        filename = self.cb.get()
        x0, y0 = self.expData.iloc[:,0], self.expData[filename]
        x_range = self.drag.getXrange()
        self.index_range = np.logical_and(x0 >= x_range[0], x0 <= x_range[1])

        self.x, self.y = x0[self.index_range], y0[self.index_range]

        try:
            if self.baseline_choose == 'auto':
                if self.alsPara is not None:
                    self.alsPara.pack_forget()
                baseline_values = peakutils.baseline(self.y)
            elif self.baseline_choose == 'ALS':
                baseline_values = self.baseline_als2(self.y,lam = float(self.lam_l.cget('text')), p = float(self.p_l.cget('text'))) if None not in (self.lam_l, self.p_l) else self.baseline_als2(self.y)

            self.line_baseline = self.ax.plot(self.x, baseline_values, label = 'baseline', color = 'red')
            self.line_corrected = self.ax.plot(self.x, self.y - baseline_values, label = 'baseline corrected', color = 'blue')
            self.ax.legend(fontsize=8).set_draggable(True)
            self.canvas.draw()
            self.menubar.entryconfig("Export", state = 'normal')
        except:
            pass


    def export_xy_files(self):
        directory = askdirectory()
        if directory !=  '':
            data_corr = self.batch()
            x = data_corr.iloc[:, 0]
            # print(data_corr)
            for i in range(len(data_corr.columns) -1):
                j = i+1
                y = data_corr.iloc[:, j]
                df = pd.DataFrame()
                df.insert(0,'x', x)
                df.insert(1, 'y', y)
                df.to_csv(os.path.join(directory, f'{self.expData.columns[j]}_corr.xy'), sep = ' ', index = False)
            messagebox.showinfo("Export to multiple .xy files", "finished!")

    def export_csv(self):
        path = asksaveasfilename(initialdir = "/",title = "Select file",filetypes = (("csv files","*.csv"),("all files","*.*")))
        data_corr = self.batch()
        data_corr.to_csv(path + '.csv', index = False)
        messagebox.showinfo("Export to a .csv file", "finished!")

    def normalization(self, y):
        return (y - min(y))/(max(y) - min(y))


    def batch(self):
        data_corr = pd.DataFrame()
        x0= self.expData.iloc[:,0]
        x_range = self.drag.getXrange()
        self.index_range = np.logical_and(x0 >= x_range[0], x0 <= x_range[1])

        self.x= x0[self.index_range]
        data_corr.insert(0, 'x', self.x)

        for i in range(len(self.expData.columns) -1):
            j = i+1
            y0 = self.expData.iloc[:, j]
            self.y= y0[self.index_range]
            y_baseline = peakutils.baseline(self.y) if self.baseline_choose == 'auto' else self.baseline_als2(self.y,lam = float(self.lam_l.cget('text')), p = float(self.p_l.cget('text')))
            y_corr = self.y - y_baseline if self.var1.get() == 0 else self.normalization(self.y - y_baseline)
            data_corr.insert(j, self.expData.columns[j], y_corr)
        return data_corr



    def polular_combobox(self, data):
        if len(data) >0:
            self.menubar.entryconfig("Import", state = 'disabled')
            self.cb.config(values = [col for col in data.columns[1:]])
            self.cb.bind("<<ComboboxSelected>>", self.callbackFunc)
            #initialize combobox
            self.cb.current(0)
            self.plotline(self.cb.get())
            #activate baseline correction menu
            self.menubar.entryconfig("Baseline", state = 'normal')

            self.set_dragRange(data)

    def set_dragRange(self, data):
        xRange = self.ax.get_xlim()
        self.drag = RangeDrag(master = self,  ax = self.ax, startPos_left = xRange[0],startPos_right = xRange[1], left_e =self.xLeft_e, right_e = self.xRight_e )

        self.canvas.draw()


    def callbackFunc(self, event):
        filename = event.widget.get()
        self.plotline(filename)
        self.baseline(self.baseline_choose)
        self.drag.update_ylim()



    def plotline(self, filename):
        # self.ax.clear()

        if self.line is not None:
            self.line.pop(0).remove()

        x, y = self.expData.iloc[:,0], self.expData[filename]
        self.line = self.ax.plot(x, y, label = f'{filename}__original', color = 'black')
        self.ax.set_ylim([min(y)-abs(max(y)*0.05), max(y)*1.05])
        self.canvas.draw()

    def baseline_als2(self, y, lam = 1e8, p =0.0001, niter=10):
        if self.alsPara is None:
            self.alsPara = LabelFrame(self, text = 'set parameters', fg = 'blue')
            self.alsPara.pack(side = 'right', anchor = 'n')
            lamV = IntVar()
            pV = IntVar()
            lam_f = LabelFrame(self.alsPara, text = 'lambda: (1e2 - 1e9)')
            p_f = LabelFrame(self.alsPara, text = 'p: (0.0001 - 0.1)')
            lam_f.pack()
            self.lam_l = Label(lam_f)
            self.p_l = Label(p_f)

            lam_s = Scale(lam_f, from_=2, to=11, resolution = 1,variable = lamV, length=600, showvalue =0, command=self.setValue_lam)
            lam_s.set(8)
            self.lam_l.config(text = 9)
            lam_f.pack(side = 'left')
            self.lam_l .pack()
            lam_s.pack()

            p_f.pack(side = 'right')
            p_s = Scale(p_f, from_=0.0001, to=0.1, resolution = 0.0001,  variable = pV, length=600,showvalue =0,  command=self.setValue_p)
            p_s.set(0.0001)
            self.p_l.config(text = 0.0001)

            p_f.pack()
            self.p_l .pack()
            p_s.pack()
        else:
            self.alsPara.pack(side = 'right', anchor = 'n')


        L = len(y)
        D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
        D = lam * D.dot(D.transpose()) # Precompute this term since it does not depend on `w`
        w = np.ones(L)
        W = sparse.spdiags(w, 0, L, L)
        for i in range(niter):
            W.setdiag(w) # Do not create a new matrix, just update diagonal values
            Z = W + D
            z = spsolve(Z, w*y)
            w = p * (y > z) + (1-p) * (y < z)
        return z

    def setValue_lam(self, val):
        self.lam_l.config(text = f'1e{val}')
        self.baseline('ALS')#update line


    def setValue_p(self, val):
        self.p_l.config(text = f'{val}')
        # print(float(self.p_l.cget('text')))
        self.baseline('ALS')#update line




def main():

    # with open('winn32') as fp:
    #     lines = fp.readlines()
    #     for line in lines:
    #         if 'qixian' in line:
    #             return


    # with open('winn32', 'r+') as fp:
    #     lines = fp.readlines()
    #     for line in lines:
    #         if '..' in line:
    #             # print(line.strip())
    #             if datetime.today().date()> datetime.strptime(line.strip().replace('..',''), '%y.%m.%d').date():
    #                 fp.write('qixian')
    #                 return

    root = Tk()
    root.title('Baseline correction')
    app = Baseline_subtraction(root)
    app.pack(fill = 'both', expand= True)
    root.mainloop()

if __name__ =='__main__':
    main()

