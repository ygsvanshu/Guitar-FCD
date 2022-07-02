import numpy as np
from PIL import Image
import matplotlib.pyplot as pl

class DataSelector:
    
    def __init__(self, ax ,l1, l2, p0, dt):
        self.ax = ax
        self.l1 = l1
        self.l2 = l2
        self.p0 = p0
        self.dt = dt
        self.x1 = list(l1.get_xdata())
        self.y1 = list(l1.get_ydata())
        self.x2 = list(l2.get_xdata())
        self.y2 = list(l2.get_ydata())
        self.x0 = list(p0.get_xdata())
        self.y0 = list(p0.get_ydata())
        self.c0 = self.ax.figure.canvas.mpl_connect('key_press_event', self.doAction)
        self.cc = self.ax.figure.canvas.mpl_connect('button_press_event', self.selectP0)
        
    def doAction(self, event):
        # print(event.key)
        if (event.key=='1'):
            self.ax.figure.canvas.mpl_disconnect(self.cc)
            self.cc = self.ax.figure.canvas.mpl_connect('button_press_event', self.selectP0)
            self.ax.set_title('Mode1')
        elif (event.key=='2'):
            self.ax.figure.canvas.mpl_disconnect(self.cc)
            self.cc = self.ax.figure.canvas.mpl_connect('button_press_event', self.selectL1)
            self.ax.set_title('Mode2')
        elif (event.key=='3'):
            self.ax.figure.canvas.mpl_disconnect(self.cc)
            self.cc = self.ax.figure.canvas.mpl_connect('button_press_event', self.selectL2)
            self.ax.set_title('Mode3')
        elif (event.key=='a'):
            alpha = np.arctan2(abs(self.x1[1]-self.x1[0]),abs(self.y1[1]-self.y1[0]))
            slope = (self.y1[1]-self.y1[0])/(self.x1[1]-self.x1[0])
            slope = -1/slope
            cross = self.y0[0] - slope*(self.x0[0])
            distn = (((self.x0[0])**2) + ((self.y0[0]-cross)**2))**0.5
            theta = np.arctan2(abs(self.x2[1]-self.x2[0]),abs(self.y2[1]-self.y2[0]))   ## FOR GENERAL USE
            # theta = np.arctan2(abs(self.y2[1]-self.y2[0]),abs(self.x2[1]-self.x2[0])) ## FOR USE IN DAY5 WITH CAMERA LINE PARALLEL INSTEAD OF PERPENDICULAR 
            self.dt += [[self.x0[0],self.y0[0],self.x1[0],self.y1[0],self.x1[1],self.y1[1],self.x2[0],self.y2[0],self.x2[1],self.y2[1],alpha,distn,theta]]
            print('n = {:2d} | α = {:6.3f} degrees | h = {:7.3f} mm | θ = {:6.3f} degrees'.format(len(self.dt),alpha*180/np.pi,distn,theta*180/np.pi))
        elif (event.key=='d'):
            try:
                removed = self.dt.pop(-1)
                alpha_r = removed[-3]
                distn_r = removed[-2]
                theta_r = removed[-1]
                print('Deleted from data: n = {:2d} | α = {:7.3f} degrees | h = {:6.3f} mm | θ = {:6.3f} degrees'.format(len(self.dt)+1,alpha_r*180/np.pi,distn_r,theta_r*180/np.pi))
            except Exception as e:
                print('Cannot delete from empty data!')
            
        self.ax.figure.canvas.draw_idle()

    def selectP0(self,event):

        if event.inaxes!=self.ax: return
        # print('click', event.xdata, event.ydata)
        if ((len(self.x0)<1)or(len(self.y0)<1)):
            self.x0.append(event.xdata)
            self.y0.append(event.ydata)
        else:
            self.x0 = []
            self.y0 = []
        self.p0.set_data(self.x0,self.y0)
        self.p0.figure.canvas.draw_idle()

    def selectL1(self, event):
        
        if event.inaxes!=self.ax: return
        if ((len(self.x1)<2)or(len(self.y1)<2)):
            self.x1.append(event.xdata)
            self.y1.append(event.ydata)
        else:
            self.x1 = []
            self.y1 = []
        self.l1.set_data(self.x1,self.y1)
        self.l1.figure.canvas.draw_idle()

    def selectL2(self, event):
        
        if event.inaxes!=self.ax: return
        if ((len(self.x2)<2)or(len(self.y2)<2)):
            self.x2.append(event.xdata)
            self.y2.append(event.ydata)
        else:
            self.x2 = []
            self.y2 = []
        self.l2.set_data(self.x2,self.y2)
        self.l2.figure.canvas.draw_idle()

class Cursor:
    """
    A cross hair cursor.
    """
    def __init__(self, ax):
        self.ax = ax
        self.horizontal_line = ax.axhline(color='k', lw=0.5, ls='--')
        self.vertical_line = ax.axvline(color='k', lw=0.5, ls='--')

    def set_cross_hair_visible(self, visible):
        need_redraw = self.horizontal_line.get_visible() != visible
        self.horizontal_line.set_visible(visible)
        self.vertical_line.set_visible(visible)
        return need_redraw

    def on_mouse_move(self, event):
        if not event.inaxes:
            need_redraw = self.set_cross_hair_visible(False)
            if need_redraw:
                self.ax.figure.canvas.draw_idle()
        else:
            self.set_cross_hair_visible(True)
            x, y = event.xdata, event.ydata
            # update the line positions
            self.horizontal_line.set_ydata(y)
            self.vertical_line.set_xdata(x)
            self.ax.figure.canvas.draw_idle()

class BlittedCursor:
    """
    A cross hair cursor using blitting for faster redraw.
    """
    def __init__(self, ax):
        self.ax = ax
        self.background = None
        self.horizontal_line = ax.axhline(color='k', lw=0.5, ls='--')
        self.vertical_line = ax.axvline(color='k', lw=0.5, ls='--')
        # text location in axes coordinates
        # self.text = ax.text(0.72, 0.9, '', transform=ax.transAxes)
        self._creating_background = False
        ax.figure.canvas.mpl_connect('draw_event', self.on_draw)

    def on_draw(self, event):
        self.create_new_background()

    def set_cross_hair_visible(self, visible):
        need_redraw = self.horizontal_line.get_visible() != visible
        self.horizontal_line.set_visible(visible)
        self.vertical_line.set_visible(visible)
        # self.text.set_visible(visible)
        return need_redraw

    def create_new_background(self):
        if self._creating_background:
            # discard calls triggered from within this function
            return
        self._creating_background = True
        self.set_cross_hair_visible(False)
        self.ax.figure.canvas.draw_idle()
        self.background = self.ax.figure.canvas.copy_from_bbox(self.ax.bbox)
        self.set_cross_hair_visible(True)
        self._creating_background = False

    def on_mouse_move(self, event):
        if self.background is None:
            self.create_new_background()
        if not event.inaxes:
            need_redraw = self.set_cross_hair_visible(False)
            if need_redraw:
                self.ax.figure.canvas.restore_region(self.background)
                self.ax.figure.canvas.blit(self.ax.bbox)
        else:
            self.set_cross_hair_visible(True)
            # update the line positions
            x, y = event.xdata, event.ydata
            self.horizontal_line.set_ydata(y)
            self.vertical_line.set_xdata(x)
            # self.text.set_text('x=%1.2f, y=%1.2f' % (x, y))

            self.ax.figure.canvas.restore_region(self.background)
            self.ax.draw_artist(self.horizontal_line)
            self.ax.draw_artist(self.vertical_line)
            # self.ax.draw_artist(self.text)
            self.ax.figure.canvas.blit(self.ax.bbox)