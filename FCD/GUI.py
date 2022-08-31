try:
    from functions import *
except:
    PRINT_ERROR('Could not import functions.py, check that the file exists in the same directory')
    exit()

try:
    import threading as td
    import tkinter as tk
    from tkinter import ttk
    from tkinter import filedialog as fd
    from tkinter import messagebox as mb
except:
    PRINT_ERROR('Could not import tkinter, try the CUI version instead of the GUI version')
    exit()

def GET_REFERENCE_IMAGE(textvar):
    textvar.set(fd.askopenfilename(title='Select reference image',initialdir='./',filetypes=[('Image files','.bmp .dds .eps .ico .jpg .jpeg .png .tga .tif .tiff')]))
    return()
    
def GET_IMAGE_DIRECTORY(textvar):
    textvar.set(fd.askdirectory(title='Select directory of images',initialdir='./'))
    return()

def SET_DATA_FILEPATH(textvar):
    tempvar = fd.asksaveasfilename(title='Save HDF5 result file',initialdir='./',filetypes=[('HDF5 files','.h5')])
    if tempvar: 
        if (tempvar[-3:]=='.h5'):
            textvar.set(tempvar)
        else:
            textvar.set(tempvar+'.h5')
    return()

def LOAD_INPUTS(inputs):
    
    filename = fd.askopenfilename(title='Save inputs as .npz file',initialdir='./',filetypes=[('Numpy zipped archive','.npz')])
    
    if (filename):
        filedata = np.load(filename)

        inputs[0].set(filedata['refimgpath'])
        inputs[1].set(filedata['imgdirpath'])
        inputs[2].set(filedata['datdirpath'])
        inputs[3].set(filedata['D'])
        inputs[4].set(filedata['units_D'])
        inputs[5].set(filedata['F'])
        inputs[6].set(filedata['units_F'])
        inputs[7].set(filedata['theta'])
        inputs[8].set(filedata['units_theta'])
        inputs[9].set(filedata['phi'])
        inputs[10].set(filedata['units_phi'])
        inputs[11].set(filedata['patternx'])
        inputs[12].set(filedata['units_patternx'])
        inputs[13].set(filedata['patterny'])
        inputs[14].set(filedata['units_patterny'])

        filename = filename.split('.')[0]
        mb.showinfo(title='Inputs loaded',message='Loaded inputs from \n{}.npz'.format(filename))

    return()

def SAVE_INPUTS(inputs):

    filename = fd.asksaveasfilename(title='Save inputs as .npz file',initialdir='./',filetypes=[('Numpy zipped archive','.npz')])

    if (filename):
        np.savez(
            filename,
            refimgpath     = inputs[0].get(),
            imgdirpath     = inputs[1].get(),
            datdirpath     = inputs[2].get(),
            D              = inputs[3].get(),
            units_D        = inputs[4].get(),
            F              = inputs[5].get(),
            units_F        = inputs[6].get(),
            theta          = inputs[7].get(),
            units_theta    = inputs[8].get(),
            phi            = inputs[9].get(),
            units_phi      = inputs[10].get(),
            patternx       = inputs[11].get(),
            units_patternx = inputs[12].get(),
            patterny       = inputs[13].get(),
            units_patterny = inputs[14].get(),
            )

        filename = filename.split('.')[0]
        mb.showinfo(title='Inputs saved',message='Current inputs have been saved to \n{}.npz'.format(filename))

    return()

def RUN_CALCULATION(inputs,progressbar,entries,comboboxes,buttons,runbutton,abortbutton,abortflag):

    def CALCERROR(message):
        mb.showerror("Error", message)
        for button in buttons:
            button.config(state='normal')
        for entry in entries:
            entry.config(state='normal')
        for combobox in comboboxes:
            combobox.config(state='readonly')
        runbutton.config(state='normal',text='Run calculation')
        abortbutton.config(state='disabled',text='! Abort !')
        progressbar.config(value=0)
        return()
    
    def CALCULATE():

        t0 = time.time()

        for button in buttons:
            button.config(state='disabled')
        for entry in entries:
            entry.config(state='disabled')
        for combobox in comboboxes:
            combobox.config(state='disabled')
        runbutton.config(state='disabled')
        abortbutton.config(state='normal')
        abortflag.set(0)

        refimgpath     = inputs[0].get()
        imgdirpath     = inputs[1].get()
        datdirpath     = inputs[2].get()
        D              = inputs[3].get()
        units_D        = inputs[4].get()
        F              = inputs[5].get()
        units_F        = inputs[6].get()
        theta          = inputs[7].get()
        units_theta    = inputs[8].get()
        phi            = inputs[9].get()
        units_phi      = inputs[10].get()
        patternx       = inputs[11].get()
        units_patternx = inputs[12].get()
        patterny       = inputs[13].get()
        units_patterny = inputs[14].get()

        linear_units  = ["m","cm","mm","ft","in"]
        linear_multi  = [1,0.01,0.001,0.3048,0.0254]
        angular_units = ["deg","rad"]
        angular_multi = [np.pi/180,1]

        if(not refimgpath):
            CALCERROR("Reference image not selected")
            return()
        if(not os.path.isfile(refimgpath)):
            CALCERROR("Reference image path is invalid")
            return()        
        if(not imgdirpath):
            CALCERROR("Directory of images not selected")
            return()        
        if(not os.path.isdir(imgdirpath)):
            CALCERROR("Directory of images is invalid")
            return() 
        if(not os.path.isdir(os.path.dirname(datdirpath))):
            CALCERROR("Path for the resulting HDF5 data file is invalid")
            return()       
        try:
            D = float(D)
        except:
            CALCERROR("D is not a number")
            return()
        if (units_D not in linear_units):
            CALCERROR("Invalid units for H")
            return()
        else:
            D = D*linear_multi[linear_units.index(units_D)]
        try:
            F = float(F)
        except:
            CALCERROR("F is not a number")
            return()
        if (units_F not in linear_units):
            CALCERROR("Invalid units for F")
            return()
        else:
            F = F*linear_multi[linear_units.index(units_F)]    
        try:
            theta = float(theta)
        except:
            CALCERROR("θ is not a number")
            return()
        if (units_theta not in angular_units):
            CALCERROR("Invalid units for θ")
            return()
        else:
            theta = theta*angular_multi[angular_units.index(units_theta)]
        try:
            phi = float(phi)
        except:
            CALCERROR("α is not a number")
            return()
        if (units_phi not in angular_units):
            CALCERROR("Invalid units for α")
            return()
        else:
            phi = phi*angular_multi[angular_units.index(units_phi)]
        try:
            patternx = float(patternx)
        except:
            CALCERROR("Pattern Δx is not a number")
            return()
        if (units_patternx not in linear_units):
            CALCERROR("Invalid units for Pattern Δx")
            return()
        else:
            patternx = patternx*linear_multi[linear_units.index(units_patternx)]
        try:
            patterny = float(patterny)
        except:
            CALCERROR("Pattern Δy is not a number")
            return()
        if (units_patterny not in linear_units):
            CALCERROR("Invalid units for Pattern Δy")
            return()
        else:
            patterny = patterny*linear_multi[linear_units.index(units_patterny)]
        imagelist  = sorted(list(filter(IS_IMAGE,os.listdir(imgdirpath))))
        if not imagelist:
            CALCERROR("No valid images in the selected directory")
            return()
        
        refimgname = refimgpath.split('/')[-1]
        refdirpath = refimgpath.replace(refimgname,'')
        if (os.path.samefile(imgdirpath,refdirpath)and(refimgname in imagelist)):
            imagelist.remove(refimgname)
        
        try:
            resultfile = h5py.File(datdirpath,'r+')
        except:
            try:
                resultfile = h5py.File(datdirpath,'w')
            except:
                CALCERROR("Unable to save HDF5 result file")
                return()
        try:
            resultfile.create_dataset('D (m)',data=[H])
        except:
            resultfile['D (m)'][0] = H
        try:
            resultfile.create_dataset('F (m)',data=[F])
        except:
            resultfile['F (m)'][0] = F
        try:
            resultfile.create_dataset('theta (rad)',data=[theta])
        except:
            resultfile['theta (rad)'][0] = theta
        try:
            resultfile.create_dataset('phi (rad)',data=[phi])
        except:
            resultfile['phi (rad)'][0] = phi

        H = D*np.cos(phi)
        c1,c2,k1,k2,nx,ny = PROCESS_IMAGE(refimgpath)
        scalex,scaley = GET_SCALE(H,F,theta,phi,k1,k2,patternx,patterny)
        SS,AT = INTGRAD2_A(nx,ny,dx=scalex,dy=scaley)

        t1 = time.time()

        print('\n[{}] Initialization complete'.format(TIMEC(t1-t0)))

        for num,image in enumerate(imagelist):

            t2 = time.time()

            imagepath = os.path.join(imgdirpath,image)
            m1,m2,tp,tp,tp,tp = PROCESS_IMAGE(imagepath)
            u,v = GET_DISPLACEMENT(m1,m2,c1,c2,k1,k2)
            h = GET_HEIGHT(A,u,v,H,theta,scalex,scaley,h00=0.0)

            imagename = image.split('.')[0]
            try:
                group = resultfile.create_group(imagename)
            except:
                group = resultfile[imagename]
            try:
                group.create_dataset('u (px)',data=u)
            except:
                group['u (px)'][:] = u
            try:
                group.create_dataset('v (px)',data=v)
            except:
                group['v (px)'][:] = v
            try:
                group.create_dataset('h (m)',data=h)
            except:
                group['h (m)'][:] = h

            t3 = time.time()

            ela = t3-t0
            eta = t3-t1
            prg = (num+1)/len(imagelist)
            eta = eta*(1-prg)/prg
            prg = int(prg*100)
            progressbar.config(value=prg)
            runbutton.config(text='Progress = {:3d}%, Elapsed {}, ETA {}'.format(prg,TIMED(ela),TIMED(eta)))
            print('[{}] Processed {}, Progress = {:3d}%, Elapsed {}, ETA {}'.format(TIMEC(t3-t2),image,prg,TIMEC(ela),TIMEC(eta)))

            if (abortflag.get()==1):
                resultfile.close()
                CALCERROR('Calculation aborted')
                print('[{}] Calculation aborted\n'.format(TIMEC(ela)))
                return()
            
        resultfile.close()
        for entry in entries:
            entry.config(state='normal')
        for button in buttons:
            button.config(state='normal')
        for combobox in comboboxes:
            combobox.config(state='readonly')
        runbutton.config(state='normal',text='Run calculation')
        abortbutton.config(state='disabled')

        t4 = time.time()
        print('[{}] Calculation completed\n'.format(TIMEC(t4-t0)))
        return()

    calc_thread = td.Thread(target=CALCULATE)
    calc_thread.start()
    return()
    
def SET_ABORT_FLAG(abort_flag,abortbutton):
    abortbutton.config(text='Aborting',state='disabled')
    abort_flag.set(1)

root = tk.Tk()
root.title('FCD Guitar')
root.resizable(False, False)

frm  = ttk.Frame(root,padding=10)
frm.grid()

refimgpath     = tk.StringVar()
imgdirpath     = tk.StringVar()
datdirpath     = tk.StringVar()

D              = tk.StringVar()
units_D        = tk.StringVar()

F              = tk.StringVar()
units_F        = tk.StringVar()

theta          = tk.StringVar()
units_theta    = tk.StringVar()

phi          = tk.StringVar()
units_phi    = tk.StringVar()

patternx       = tk.StringVar()
units_patternx = tk.StringVar()

patterny       = tk.StringVar()
units_patterny = tk.StringVar()

abort_flag     = tk.IntVar()
abort_flag.set(0)

l1 = ttk.Label(frm,text="Reference image",anchor=tk.E)
l1.grid(column=0,row=0,sticky='nsew')

e1 = ttk.Entry(frm,textvariable=refimgpath)
e1.grid(column=1,row=0,columnspan=5,sticky='nsew')

b1 = ttk.Button(frm,text="Browse", command=lambda:GET_REFERENCE_IMAGE(refimgpath))
b1.grid(column=6,row=0,columnspan=2,sticky='nsew')

l2 = ttk.Label(frm,text="Directory of images",anchor=tk.E)
l2.grid(column=0,row=1,sticky='nsew')

e2 = ttk.Entry(frm,textvariable=imgdirpath)
e2.grid(column=1,row=1,columnspan=5,sticky='nsew')

b2 = ttk.Button(frm, text="Browse", command=lambda:GET_IMAGE_DIRECTORY(imgdirpath))
b2.grid(column=6,row=1,columnspan=2,sticky='nsew')

l3 = ttk.Label(frm,text="Save result HDF5 file",anchor=tk.E)
l3.grid(column=0,row=2,sticky='nsew')

e3 = ttk.Entry(frm,textvariable=datdirpath)
e3.grid(column=1,row=2,columnspan=5,sticky='nsew')

b3 = ttk.Button(frm,text="Browse", command=lambda:SET_DATA_FILEPATH(datdirpath))
b3.grid(column=6,row=2,columnspan=2,sticky='nsew')

s1 = ttk.Separator(frm,orient='horizontal')
s1.grid(column=0,row=3,columnspan=7,sticky='nsew',pady=10)

s2 = ttk.Separator(frm,orient='vertical')
s2.grid(column=3,row=4,rowspan=3,sticky='nsew',padx=10)

l4 = ttk.Label(frm,text="D = ",anchor=tk.E)
l4.grid(column=0,row=4,sticky='nsew')

e4 = ttk.Entry(frm,textvariable=D)
e4.grid(column=1,row=4,columnspan=1,sticky='nsew')

d1 = ttk.Combobox(frm,values=["m","cm","mm","ft","in"],textvariable=units_D,exportselection=0,state='readonly')
d1.grid(column=2,row=4,sticky='nsew')

l5 = ttk.Label(frm,text="F = ",anchor=tk.E)
l5.grid(column=4,row=4,sticky='nsew')

e5 = ttk.Entry(frm,textvariable=F)
e5.grid(column=5,row=4,columnspan=1,sticky='nsew')

d2 = ttk.Combobox(frm,values=["m","cm","mm","ft","in"],textvariable=units_F,exportselection=0,state='readonly')
d2.grid(column=6,row=4,sticky='nsew')

l6 = ttk.Label(frm,text="θ = ",anchor=tk.E)
l6.grid(column=0,row=5,sticky='nsew')

e6 = ttk.Entry(frm,textvariable=theta)
e6.grid(column=1,row=5,columnspan=1,sticky='nsew')

d3 = ttk.Combobox(frm,values=["deg","rad"],textvariable=units_theta,exportselection=0,state='readonly')
d3.grid(column=2,row=5,sticky='nsew')

l7 = ttk.Label(frm,text="ϕ = ",anchor=tk.E)
l7.grid(column=4,row=5,sticky='nsew')

e7 = ttk.Entry(frm,textvariable=phi)
e7.grid(column=5,row=5,columnspan=1,sticky='nsew')

d4 = ttk.Combobox(frm,values=["deg","rad"],textvariable=units_phi,exportselection=0,state='readonly')
d4.grid(column=6,row=5,sticky='nsew')

l8 = ttk.Label(frm,text="pattern Δx = ",anchor=tk.E)
l8.grid(column=0,row=6,sticky='nsew')

e8 = ttk.Entry(frm,textvariable=patternx)
e8.grid(column=1,row=6,columnspan=1,sticky='nsew')

d5 = ttk.Combobox(frm,values=["m","cm","mm","ft","in"],textvariable=units_patternx,exportselection=0,state='readonly')
d5.grid(column=2,row=6,sticky='nsew')

l9 = ttk.Label(frm,text="pattern Δy = ",anchor=tk.E)
l9.grid(column=4,row=6,sticky='nsew')

e9 = ttk.Entry(frm,textvariable=patterny)
e9.grid(column=5,row=6,columnspan=1,sticky='nsew')

d6 = ttk.Combobox(frm,values=["m","cm","mm","ft","in"],textvariable=units_patterny,exportselection=0,state='readonly')
d6.grid(column=6,row=6,sticky='nsew')

s3 = ttk.Separator(frm,orient='horizontal')
s3.grid(column=0,row=7,columnspan=7,sticky='nsew',pady=10)

s4 = ttk.Separator(frm,orient='vertical')
s4.grid(column=3,row=8,rowspan=2,sticky='nsew',padx=10)

b4 = ttk.Button(frm,text='Load inputs',state='normal')
b4.grid(column=0,row=8,columnspan=3,sticky='nsew')

b5 = ttk.Button(frm,text='Save inputs ',state='normal')
b5.grid(column=0,row=9,columnspan=3,sticky='nsew')

p1 = ttk.Progressbar(mode='determinate')
p1.grid(column=0,row=10,sticky='nsew',padx=10,pady=10)

b6 = ttk.Button(frm,text='Run calculation',state='normal',width=50)
b6.grid(column=4,row=8,columnspan=3,sticky='nsew')

b7 = ttk.Button(frm,text='! Abort !',state='disabled')
b7.grid(column=4,row=9,columnspan=3,sticky='nsew')

inputs  = [refimgpath,imgdirpath,datdirpath,D,units_D,F,units_F,theta,units_theta,phi,units_phi,patternx,units_patternx,patterny,units_patterny]
entries = [e1,e2,e3,e4,e5,e6,e7,e8,e9]
comboboxes = [d1,d2,d3,d4,d5,d6]
buttons = [b1,b2,b3,b4,b5]

b4.config(command=lambda:LOAD_INPUTS(inputs))
b5.config(command=lambda:SAVE_INPUTS(inputs))
b6.config(command=lambda:RUN_CALCULATION(inputs,p1,entries,comboboxes,buttons,b6,b7,abort_flag))
b7.config(command=lambda:SET_ABORT_FLAG(abort_flag,b7))

root.mainloop()

exit()