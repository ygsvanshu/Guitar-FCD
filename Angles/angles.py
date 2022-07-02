from functions import *
        
im = Image.open("PNGs/guitarday5.png")
sx = 297
sy = 420

print(im.format, im.size, im.mode)
im = np.array(im.convert('1'))

# data = [[220.14107195,60.51910899,237.79111941,42.54439511,220.70185941,59.84026101,197.42721694,285.93566796,295.69861743,194.17157682,0.79140667,174.36458265,0.81962739]]
data = []

fig, ax = pl.subplots()
ax.set_title('Mode1')
ax.plot([0,1,2,3],[0,1,2,3],'k-')
# ax.imshow(im,cmap='binary_r',extent=[0,sx,0,sy])
l1, = ax.plot([],[],'r:.')  # empty line
l2, = ax.plot([],[],'g:.')  # empty line
p0, = ax.plot([],[],'b.')  # empty line
linebuilder = DataSelector(ax,l1,l2,p0,data)
cursor = BlittedCursor(ax)
fig.canvas.mpl_connect('motion_notify_event', cursor.on_mouse_move)
pl.show()

filename = input('Save data file as: ')
if filename: np.save(filename,np.array(data))