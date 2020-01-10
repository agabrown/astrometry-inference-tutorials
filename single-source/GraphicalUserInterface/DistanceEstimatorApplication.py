# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 13:28:04 2017

@author: Ariadna
"""

import tkinter as tk
from tkinter import ttk, filedialog
from main import *
from pyrallaxes import *
import numpy as np
from PIL import ImageTk, Image
import pickle
import warnings
#import Image as im
import matplotlib.pyplot as plt
try:
    # NavigationToolbar2TkAgg was deprecated in matplotlib 2.2 and is no
    # longer present in v3.0 - so a version test would also be possible.
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
except ImportError:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
except:
    raise
from matplotlib.figure import Figure
#from main import *
LARGE_FONT = ("Verdana", 12)


class WarningWindow(tk.Tk):
    def __init__(self,message,*args,**kwargs):
        tk.Tk.__init__(self,*args,**kwargs)
        tk.Tk.wm_title(self,'WARNING!')
        self.container = tk.Frame(self)
        self.container.pack(side="top", fill="both", expand=True)
        self.container.grid_rowconfigure(0, weight=1)
        self.container.grid_columnconfigure(0, weight=1)
        
        lb = tk.Label(self.container,text=message)
        lb.pack()
        

class PDFWindow(tk.Tk):
    def __init__(self,f,w, s,par,r1,r2,r3,r4,r0,par2,n,st,st2,*args,**kwargs):
        tk.Tk.__init__(self,*args,**kwargs)
        tk.Tk.wm_title(self,st2)
        
        self.container = tk.Frame(self)
        self.container.pack(side="top", fill="both", expand=True)
        self.container.grid_rowconfigure(0, weight=1)
        self.container.grid_columnconfigure(0, weight=1)
        
        
        
        self.plot_image(f,w, s,par,r1,r2,r3,r4,r0,par2,n,st)
        
        
        button = tk.Button(self.container,text = 'Save Image',command = lambda: self.save_image())
        button.pack(padx=100,pady=10)
        
        self.container.tkraise()
        
        
        
        
        
    
    
        
    def plot_image(self,f,w, s,par,r1,r2,r3,r4,r0,par2,n,st):
        #print('in plot_image',r3,r1,r2,r4,n)
        self.fig = Figure(figsize=(5,4), dpi=100)
        pdf = self.fig.add_subplot(111)
        
        t = np.arange(r0,par2,par2/500)
        y = [f(x,w,s,par)/n for x in t]
        
        pdf.plot(t,y)
        
        plt.show()
        dim = f(r1,w,s,par)/n
        t2 = r1*np.ones(20)
        y2 = np.arange(0,dim,dim/20)
        pdf.plot(t2,y2,label = 'Mode',color = 'red')
        
        dim_med = f(r2,w,s,par)/n
       
        if dim_med != 0:
            t3 = r2*np.ones(20)
            y3 = np.arange(0,dim_med,dim_med/20)
        
            pdf.plot(t3,y3,label='Median',color = 'black')
        
        
        x = np.arange(r3,r4,1)
        p = [f(i,w,s,par)/n for i in x]
        
        pdf.fill_between(x,0,p,facecolor = 'blue',alpha = 0.5,label = '90 % uncertainty interval')
        
        pdf.legend()
        if st == 'r':
            pdf.set_xlabel('True distance r (kpc)')
        else:
            pdf.set_xlabel('Distance Modulus')
        pdf.set_ylabel('PDF')
        pdf.set_title('Probability Distribution Function')
        
        dataPlot = FigureCanvasTkAgg(self.fig, master=self.container)
        if hasattr(dataPlot, 'show'):
            # Deprecated version
            dataPlot.show()
        elif hasattr(dataPlot, 'draw'):
            dataPlot.draw()
        else:
            warnings.warn("Cannot show the plot with either dataPlot.show() or dataPlot.draw()")
            return
        dataPlot.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        
        
        
    
    def save_image(self):
        image = self.fig
        ftypes = [('Image files', '*.png'), ('All files', '*')]
        dlg = filedialog.SaveAs(self,defaultextension='.png',filetypes = ftypes)
        fl = dlg.show()
        
        if fl != '':
            if not fl.endswith('.png'):
                fl = fl +'.png'                        
            image.savefig(fl)
        

class DistanceEstimator(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        tk.Tk.wm_title(self, "Distance Estimation Tool")
        
        self.container = tk.Frame(self)
        self.container.pack(side="top", fill="both", expand=True)
        self.container.grid_rowconfigure(0, weight=1)
        self.container.grid_columnconfigure(0, weight=1)

        self.frames = {}

        for F in (StartPage, PageOne, PageTwo, PageThree):
            frame = F(self.container, self)

            self.frames[F] = frame

            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(StartPage)
        
    
    def compute_results(self,cont):
        r_lim = np.float64(100) #kpc
        r0 = 0.01 #kpc
        L = 1.35 #kpc
        m0 = -10 
        mlim = 5*np.log(100*1000)-5        
       
        frame = self.frames[cont]
        
        w = float(self.frames[PageOne].parallax.get())
        s = float(self.frames[PageOne].parallax_err.get())
        
        
        # Defining the booleans corresponding to the options marked by the users i the checkbuttons
        boolun = self.frames[PageOne].un.get()
        boolexp = self.frames[PageOne].exp.get()
        booltrans = self.frames[PageOne].trans.get()
        
        boolmun = self.frames[PageOne].mun.get()
        boolmexp = self.frames[PageOne].mexp.get()
        boolmtrans = self.frames[PageOne].mtrans.get()
        
        # Computing the mode, the median, the 5 % and the 95 % percentiles by the selected method 
        
        if boolun == True:
            r_ud_5,r_mode_ud,r_media_ud,r_ud_95,n_ud = main_un(w,s)
            
        if boolexp == True:
            r_exp_5,r_mode_exp,r_media_exp,r_exp_95,n_exp = main_exp(w,s)
            
        if booltrans == True:
            rsmin,rs,rsmax = main_trans(w,s)
            
        if boolmun == True:
            r_mud_5,m_mode_un,r_media_mud,r_mud_95,n_mud = main_mun(w,s)
            
        if boolmexp == True:
            r_mexp_5,m_exp,r_media_mexp,r_mexp_95,n_mexp = main_mexp(w,s) 
            
        if boolmtrans == True:
            mtm_min,mtm,mtm_max = main_m_trans(w,s)
           
        
        
        
        lf1 = tk.LabelFrame(frame,text = 'Legend')
        lf1.pack(padx = 5, pady = 10)
        
        index = 0
        if boolun == True:
            index = index + 1
            
            lab1 = tk.Label(lf1,text = 'r_UDP = ')
            lab1.grid(row = index, column = 0)
            
            lab11 = tk.Label(lf1,text = 'Distance using Uniform Distance Prior')
            lab11.grid(row = index, column = 1)
 
        if boolexp == True:
            index = index + 1
            
            lab2 = tk.Label(lf1,text = 'r_EDSDP = ')
            lab2.grid(row = index, column = 0)
            
            lab22 = tk.Label(lf1,text = 'Distance using Exponentially Decreasing Space Density Prior ')
            lab22.grid(row = index, column = 1)

        
        
        if boolmun == True:
            index = index + 1
            
            lab4 = tk.Label(lf1,text = 'm_UDP = ')
            lab4.grid(row = index, column = 0)
            
            lab44 = tk.Label(lf1,text = 'Distance Modulud using Uniform Distance Prior')
            lab44.grid(row = index, column = 1)
            
        if boolmexp == True:
            index = index + 1
            
            lab5 = tk.Label(lf1,text = 'm_EDSDP = ')
            lab5.grid(row = index, column = 0)
            
            lab55 = tk.Label(lf1,text = 'Distance Modulus using Exponentially Decreasing Space Density Prior')
            lab55.grid(row = index, column = 1)
            
        if booltrans == True:            
            index = index + 1
            
            lab3 = tk.Label(lf1,text = 'r_TM = ')
            lab3.grid(row = index, column = 0)
             
            lab33 = tk.Label(lf1,text = 'Distance using Transformation Method')
            lab33.grid(row = index, column = 1)
            
        if boolmtrans == True:
            index = index + 1
            
            lab6 = tk.Label(lf1,text = 'm_TM = ')
            lab6.grid(row = index, column = 0)
            
            lab66 = tk.Label(lf1,text = ' Distance Modulus using Transformation Method')
            lab66.grid(row = index, column = 1)
            
        if (boolun or boolexp or boolmun or boolmexp) == True: 
            
            lf2 = tk.LabelFrame(frame,text = 'Bayesian Methods')
            lf2.pack(padx = 5, pady = 10)
            
            l1 = tk.Label(lf2,text = 'Mode (pc)')
            l1.grid(row = 0,column = 1)
            
            l2 = tk.Label(lf2,text = 'Median (pc)')
            l2.grid(row = 0,column = 2)
            
            l3 = tk.Label(lf2,text = '5% quantile (pc)')
            l3.grid(row = 0,column = 3)
            
            l4 = tk.Label(lf2,text = '95% quantile (pc)')
            l4.grid(row = 0,column = 4)      
            
        lf3 = tk.LabelFrame(frame,text ='Transform Method')
        lf3.pack(padx = 5, pady = 10)
            
        if (booltrans or boolmtrans) == True:
            
            l5 = tk.Label(lf3,text = 'Estimate (pc)')
            l5.grid(row = 0, column = 1)
            
            l6 = tk.Label(lf3,text = 'Inferior bound (pc)')
            l6.grid(row = 0, column = 2)
            
            l7 = tk.Label(lf3,text = 'Superior bound (pc)')
            l7.grid(row = 0, column = 3)
            
        index = 0
            
        if boolun == True:
            index = index + 1
            l10 = tk.Label(lf2,text = 'r_UDP')
            l10.grid(row = index,column = 0)
                
            l11 = tk.Label(lf2,text = str(r_mode_ud))
            l11.grid(row = index,column = 1)
            
            l12 = tk.Label(lf2,text = str(r_media_ud))
            l12.grid(row = index,column = 2)
            
            l13 = tk.Label(lf2,text = str(r_ud_5))
            l13.grid(row = index,column = 3)
            
            l14 = tk.Label(lf2,text = str(r_ud_95))
            l14.grid(row = index,column = 4)            
            
            button1 = tk.Button(lf2,text = 'Show PDF', command = lambda: self.show_figure(uniform_distance_posterior,w,s,r_lim,r_mode_ud/1000,r_media_ud/1000,r_ud_5/1000,r_ud_95/1000,r0,r_lim,n_ud,'r','Uniform Distance PDF'))            
            button1.grid(row = index,column = 5)
            
        if boolexp == True:
            index = index + 1
            l20 = tk.Label(lf2,text = 'r_EDSDP')
            l20.grid(row = index,column = 0)
                
            l21 = tk.Label(lf2,text = str(r_mode_exp))
            l21.grid(row = index,column = 1)
            
            
            l22 = tk.Label(lf2,text = str(r_media_exp))
            l22.grid(row = index,column = 2)
            
            l23 = tk.Label(lf2,text = str(r_exp_5))
            l23.grid(row = index,column = 3)
            
            l24 = tk.Label(lf2,text = str(r_exp_95))
            l24.grid(row = index,column = 4)            
            
            button2 = tk.Button(lf2,text = 'Show PDF', command = lambda: self.show_figure(exponentially_decreasing_space_density_posterior,w,s,L,r_mode_exp/1000,r_media_exp/1000,r_exp_5/1000,r_exp_95/1000,r0,r_lim,n_exp,'r','Exponentially Decreasing Space Density PDF'))            
            button2.grid(row = index,column = 5)
            
        if boolmun == True:
            index = index + 1
            l30 = tk.Label(lf2,text = 'm_UDP')
            l30.grid(row = index,column = 0)
                
            l31 = tk.Label(lf2,text = str(m_mode_un))
            l31.grid(row = index,column = 1)
                
            l32 = tk.Label(lf2,text = str(r_media_mud))
            l32.grid(row = index,column = 2)
            
            l33 = tk.Label(lf2,text = str(r_mud_5))
            l33.grid(row = index,column = 3)
            
            
            l34 = tk.Label(lf2,text = str(r_mud_95))
            l34.grid(row = index,column = 4)            
            
            button3 = tk.Button(lf2,text = 'Show PDF', command = lambda: self.show_figure(dmpdfun,w/1000,s/1000,r_lim*1000,m_mode_un,r_media_mud,r_mud_5,r_mud_95,m0,mlim,n_mud,'m','Uniform Distance Modulus PDF'))            
            button3.grid(row = index,column = 5)
            
            
        if boolmexp == True:
            index = index + 1
            l40 = tk.Label(lf2,text = 'm_EDSDP')
            l40.grid(row = index,column = 0)
                
        
            l41 = tk.Label(lf2,text = str(m_exp))
            l41.grid(row = index,column = 1)
                
            l42 = tk.Label(lf2,text = str(r_media_mexp))
            l42.grid(row = index,column = 2)
            
            l43 = tk.Label(lf2,text = str(r_mexp_5))
            l43.grid(row = index,column = 3)
            
            l44 = tk.Label(lf2,text = str(r_mexp_95))
            l44.grid(row = index,column = 4)            
            
            button4 = tk.Button(lf2,text = 'Show PDF', command = lambda: self.show_figure(dmpdfexp,w/1000,s/1000,L*1000,m_exp,r_media_mexp,r_mexp_5,r_mexp_95,m0,mlim,n_mexp,'m','Exponentially Decreasing Space Density Distance Modulus PDF'))            
            button4.grid(row = index,column = 5)
            
        index = 0
        if booltrans == True:
            index = index + 1
                
            l50 = tk.Label(lf3,text = 'r_TM')
            l50.grid(row = index,column = 0)
                
            l51 = tk.Label(lf3,text = str(rs))
            l51.grid(row = index,column = 1)
                
            l52 = tk.Label(lf3,text = str(rsmin))
            l52.grid(row = index,column = 2)
            
            l53 = tk.Label(lf3,text = str(rsmax))
            l53.grid(row = index,column = 3)
                
        if boolmtrans == True:
            index = index + 1
                
            l60 = tk.Label(lf3,text = 'm_TM')
            l60.grid(row = index,column = 0)
                
            l61 = tk.Label(lf3,text = str(mtm))
            l61.grid(row = index,column = 1)
                
            l62 = tk.Label(lf3,text = str(mtm_min))
            l62.grid(row = index,column = 2)
            
            l63 = tk.Label(lf3,text = str(mtm_max))
            
            l63.grid(row = index,column = 3)
            
            
                
    def show_frame(self, cont):
        
        frame = self.frames[cont]
        frame.tkraise()
        
    def show_frame_2(self,cont):
        
        w = self.frames[PageOne].parallax.get()
        s = self.frames[PageOne].parallax_err.get()
        
        if w != "" and s != "":
            frame = cont(self.container,self)
            self.frames[cont] = frame
            frame.grid(row=0, column=0, sticky="nsew")
        
            self.compute_results(cont)
            frame.tkraise()
            
        else:
            message = 'It is mandatory to introduce an observed parallax and its parallax error in miliarcseconds!'
            warning = WarningWindow(message)
            warning.mainloop()
         
        
        
    def show_figure(self,f,w, s,par,r1,r2,r3,r4,r0,par2,n,st,st2):
        
        fig = PDFWindow(f,w, s,par,r1,r2,r3,r4,r0,par2,n,st,st2)
        fig.mainloop()
    
    
        
    
            
            
    
        
class StartPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Distance estimator", font=LARGE_FONT)
        label.pack(pady=10, padx=10)

        button = ttk.Button(self, text="Input Data",
                            command=lambda: controller.show_frame(PageOne))

        button.pack(side=tk.TOP,padx=5,pady=10)

        button2 = ttk.Button(self, text="Info",
                             command=lambda: controller.show_frame(PageTwo))
        button2.pack(side=tk.TOP)
        
        path = 'log.png'
        im = Image.open(path)
        tkimage = ImageTk.PhotoImage(im, master=self)
        img = tk.Label(self, image = tkimage)
        img.image = tkimage
        img.pack(side=tk.LEFT,padx = 10,pady=10)
        
        text = 'This tool has been developed in the framework of the collaboration scholarship AGAUR with the "Departament de Física Quàntica i Astrofísica" of the University of Barcelona. \n \n The Distance Estimator Tool offers the possibility of computing distances and distance modulus from trigonometric parallaxes using both Bayesian Methods and a frequentist method called the Transformation Method.\n \n Given an observed parallax and its associated error it provides the mode and the median of the true distance distribution as well as the 90% uncertainty interval. In the case of the distance modulus using the Transformation Method you can also provide the observed distance modulus. If no number is introduced the default is 0.'
        label1 = tk.Label(self,text = text,wraplength = 300,padx = 10,pady=5,font = 16,justify = tk.LEFT)
        label1.pack(side=tk.RIGHT)
        
       

class PageOne(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Input Data", font=LARGE_FONT)
        label.pack(pady=10, padx=10)

        button1 = ttk.Button(self, text="Back to Home",
                             command=lambda: controller.show_frame(StartPage))
        button1.pack()

        button2 = ttk.Button(self, text="Confirm Selected Options",
                             command=lambda: controller.show_frame_2(PageThree))
        button2.pack(side=tk.BOTTOM,padx=10,pady=20)
        
        
        lf1 = ttk.LabelFrame(self,text='Input data')
        lf1.pack(side=tk.TOP,fill=tk.X,padx=5,pady=30)
        
        frame1 = tk.Frame(lf1)
        frame1.pack()
        
        lbl1 = tk.Label(frame1, text="Observed parallax (mas)", width=20)
        lbl1.grid(row = 0,column = 0)
        
        self.parallax = tk.StringVar()
        self.parallax.set("")
        
        self.entry0 = tk.Entry(frame1)
        self.entry0["textvariable"] = self.parallax
        self.entry0.grid(row = 0, column = 1)

        lbl2 = tk.Label(frame1,text= 'Parallax error (mas)',width = 20)
        lbl2.grid(row = 1, column = 0)
        self.parallax_err = tk.StringVar()
        self.parallax_err.set("")
        
        entry1 = tk.Entry(frame1)
        entry1["textvariable"] =self.parallax_err
        entry1.grid(row = 1, column = 1)
       
        
        lf2 = ttk.LabelFrame(self,text='Distance')
        lf2.pack(side=tk.LEFT,fill=tk.X,padx=5,pady=30)
        
        self.un = tk.BooleanVar()
        self.exp = tk.BooleanVar()
        self.trans = tk.BooleanVar()
        
        

        cb1 = tk.Checkbutton(lf2, text="Uniform Distance Prior",
                             variable=self.un)
        cb1.select()
        cb1.place(x=10, y=50)
        cb1.pack()

        cb2 = tk.Checkbutton(lf2, text="Exponentially Decreasing Space Density Prior",
                             variable=self.exp)
        cb2.select()

        cb2.pack( padx=10, pady=10)

        cb3 = tk.Checkbutton(lf2, text="Transformation Method",
                             variable=self.trans)
        cb3.select()

        cb3.pack( padx=10, pady=10)
        
        lf3 = ttk.LabelFrame(self,text='Distance Modulus')
        lf3.pack(side=tk.RIGHT,fill=tk.X,padx=5,pady=30)
        
        self.mun = tk.BooleanVar()
        self.mexp = tk.BooleanVar()
        self.mtrans = tk.BooleanVar()
        
        
        
        cb4 = tk.Checkbutton(lf3, text="Uniform Distance Prior",
                             variable=self.mun)
        cb4.select()
        cb4.place(x=10, y=50)
        cb4.pack()

        cb5 = tk.Checkbutton(lf3, text="Exponentially Decreasing Space Density Prior",
                             variable=self.mexp)
        cb5.select()

        cb5.pack( padx=10, pady=10)

        cb6 = tk.Checkbutton(lf3, text="Transformation Method",
                             variable=self.mtrans)
        cb6.select()

        cb6.pack( padx=10, pady=10)
        
        
     


class PageTwo(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Additional information", font=LARGE_FONT)
        label.pack(pady=10, padx=10)

        button1 = ttk.Button(self, text="Back to Home",
                             command=lambda: controller.show_frame(StartPage))
        button1.pack()

        button2 = ttk.Button(self, text="Input Data",
                             command=lambda: controller.show_frame(PageOne))
        button2.pack()
        
        text = " "
        lf = tk.Label(self,text = text,wraplength = 300,padx = 10,pady=5,font = 16,justify = tk.LEFT)
        lf.pack()
        

class PageThree(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Distance and distance modulus estimates", font=LARGE_FONT)
        label.pack(pady=10, padx=10)
        self.controller = controller
        button1 = ttk.Button(self, text="Back to Home",
                             command=lambda: controller.show_frame(StartPage))
        button1.pack()
        
       
        
        
        
        
       
        

app = DistanceEstimator()
app.mainloop()
