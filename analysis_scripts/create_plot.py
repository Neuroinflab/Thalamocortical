from lfp_parameters import *
from matplotlib.mlab import griddata
from matplotlib import rcParams
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


rcParams.update({'font.size': 8})
def update_contour_plot(i, data,  ax, xi, yi, levels):
    ax.cla()
    im = ax.contourf(xi.reshape(num_x, num_y), 
                     yi.reshape(num_x, num_y),  
                     data[:,i].reshape(num_x, num_y), 15, 
                     levels = levels, cmap=plt.cm.PRGn)
    plt.title('Time='+str(i))
    plt.ylim(ymin=-2400,ymax=550)
    plt.xlim(xmin=-450,xmax=450)
    return im,

def grid(x, y, z, resX=100, resY=100):
    "Convert 3 column data to matplotlib grid"
    z = z.flatten()
    xi = np.linspace(min(x), max(x), resX)
    yi = np.linspace(min(y), max(y), resY)
    Z = griddata(x, y, z, xi, yi)
    X, Y = np.meshgrid(xi, yi)
    return X, Y, Z

lfp = np.load(pot_filename)
fig = plt.figure()#figsize=(4,6))
plt.subplot(121, aspect='equal')
for pop_name in pop_names:
    src_pos = fetch_mid_pts(h, pop_name)
plt.scatter(src_pos[:, 0], src_pos[:, 1], marker='.', alpha=1, color='k', lw = 1., s=0.5)
plt.scatter(ele_pos[:, 0], ele_pos[:, 1], marker='o', alpha=1., color='r', lw = 1., s=0.4)
plt.xlabel('X ($\mu$m)')
plt.ylabel('Y ($\mu$m)')
plt.title('Morphology and electrodes')
plt.ylim(ymin=-2400,ymax=550)
plt.xlim(xmin=-450,xmax=450)

if ele_config == '1D':
    ax = plt.subplot(122)
    # for ele in range(num_ele):
    #     plt.plot(lfp[ele, :]+ele/2.)
    #ax.yaxis.set_ticks(ele_pos[:, 1])
    plt.imshow(lfp[::-1], aspect='auto', interpolation='nearest', cmap=plt.cm.PRGn)

    plt.yticks(range(num_ele), ele_pos[:, 1].astype(int)[::-1])
    plt.xticks(np.arange(2500, 4500, 500), np.arange(250, 450, 50))
    plt.ylabel('Electrode depth ($\mu$m)')
    plt.xlabel('Time (ms)')
    plt.title('Extracellular potentails')
    #plt.xlim(xmin=2500, xmax=4500)
    plt.colorbar()

elif ele_config == '2D':
    ax2 = plt.subplot(122, aspect='equal')
    levels = np.linspace(-2.,2.,15)
    ### Contour plot at an instance of time!
    im2 = plt.contourf(ele_pos[:,0].reshape(num_x, num_y), 
                       ele_pos[:,1].reshape(num_x, num_y), 
                       lfp[:,1105].reshape(num_x,num_y), 15, 
                       levels = levels, cmap=plt.cm.PRGn)
    plt.colorbar()
    plt.title('Time=110.5ms')
    plt.ylim(ymin=-2400,ymax=550)
    plt.xlim(xmin=-450,xmax=450)
    if ANIMATE:
        ani = animation.FuncAnimation(fig, update_contour_plot, frames=xrange(1100,1200), 
                                      fargs=(lfp, ax2, ele_pos[:,0], ele_pos[:,1], levels), interval=50)
plt.tight_layout()
#plt.subplots_adjust(wspace = 0.001)
#plt.setp(ax2.xaxis.get_majorticklabels(), rotation=70 )
#plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
plt.savefig('1D_1105_traub.png', dpi=600)
#plt.show()
h.close()
