import scipy.spatial
from scipy import signal
import numpy as np
import h5py as h5
import brewer2mpl
from matplotlib import rcParams, cm, gridspec
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes

rcParams.update({'font.size': 8, 'font.family': 'sans-serif'})
def fetch_mid_pts(h, pop_name):
    all_pts = h['/data/static/morphology/'+pop_name]
    x = (all_pts['x0']+all_pts['x1']) / 2.
    y = (all_pts['y0']+all_pts['y1']) / 2.
    z = (all_pts['z0']+all_pts['z1']) / 2.
    x = x.reshape(x.size,1)
    y = y.reshape(y.size,1)
    z = z.reshape(z.size,1)
    return np.hstack((x, y, z))

def find_map_idx(h, pop_name, field_name):
    '''Find the corresponding, locations of morphology in field'''
    mor = h['/data/static/morphology/'+pop_name]
    i_data = h['/data/uniform/'+pop_name+'/'+field_name]
    mor_names = mor.dims[0].values()[0] #entry in /map
    i_names = i_data.dims[0].values()[0] #entry in /map
    if np.array_equal(i_names, mor_names):
        idx = range(i_names.shape[0])
    else:
        idx = [np.where(i_names.value==entry)[0][0] for entry in mor_names]
    return idx

def inv_distance(src_pos, ele_pos):
    '''computes the inverse distance between src_pos and ele_pos'''
    dist_matrix = np.zeros((src_pos.shape[0], ele_pos.shape[0]))
    for ii,electrode in enumerate(ele_pos):
        dist_matrix[:, ii] = scipy.spatial.distance.cdist(src_pos, electrode.reshape(1,3)).flatten()
    dist_matrix = 1 / dist_matrix #inverse distance matrix
    return dist_matrix

def pot_vs_time(h, pop_name, field_name, src_pos, ele_pos):
    '''returns potentials, at ele_pos, due to src_pos, over time'''
    idx = find_map_idx(h, pop_name, field_name)
    src_time = h['/data/uniform/'+pop_name+'/'+field_name].value
    src_time = src_time[idx] # Order according to correct indices
    ele_src = inv_distance(src_pos, ele_pos).T
    return np.dot(ele_src, src_time)*(1 / (4*np.pi*0.3))

def get_all_src_pos(h, pop_names, total_cmpts):
    all_srcs = np.zeros((sum(total_cmpts), 3))
    for jj, pop_name in enumerate(pop_names):
        all_srcs[np.sum(total_cmpts[:jj]):np.sum(total_cmpts[:jj+1]), :] = fetch_mid_pts(h, pop_name)
    return all_srcs

def get_extracellular(h, pop_names, time_pts, ele_pos, variable):
    pot_sum = np.zeros((num_ele, time_pts))
    for pop_name in pop_names:
        src_pos = fetch_mid_pts(h, pop_name)
        pot_sum += pot_vs_time(h, pop_name, variable, src_pos, ele_pos)
        print 'Done extracellular pots for pop_name, using currents', pop_name, variable
    return pot_sum

def place_electrodes_1D(n):
    '''places n number of electrodes in ID along the column'''
    ele_x = np.ones((n,1))*25.
    ele_y = np.linspace(-2050., 450.0, num=n).reshape(n,1)
    ele_z = np.ones((n,1))*25.
    return np.hstack((ele_x, ele_y, ele_z))

def lowpassfilter(data):
    sampling_rate = 10000 #10kHz, corresp to 0.1ms
    nyq = 0.5 * sampling_rate
    normal_cutoff = 100 / nyq
    b, a = signal.butter(2, normal_cutoff, 'low', analog=False)
    return signal.lfilter(b, a, data)

def plot_potential(ax, lfp, max_val, title):
    norm = cm.colors.Normalize(vmax=max_val, vmin=-max_val, clip=False)
    im = plt.imshow(lfp[::-1], aspect='auto', norm=norm, interpolation='nearest', cmap=plt.cm.PRGn)
    plt.xlim((2750, 4250))
    plt.xticks(np.arange(3000, 5000, 1000), np.arange(300, 500, 100))
    #plt.ylabel('Electrode depth ($\mu$m)')
    #plt.xlabel('Time (ms)')
    plt.title(title, fontweight="bold", fontsize=12)
    #plt.xlim(xmin=2500, xmax=4500)
    #plt.colorbar(extend='both')
    cbaxes = inset_axes(ax,
                        width="40%",  # width = 10% of parent_bbox width
                        height="3%",  # height : 50%
                        loc=1, borderpad=1)
    cbar = plt.colorbar(cax=cbaxes, ticks=[-max_val,0.,max_val], orientation='horizontal', format='%.2f')
    cbar.ax.set_xticklabels([round(-max_val,2),0.,round(max_val,2)])
    return ax, im

def get_specific(ax, h, pop_names, time_pts, ele_pos, variable, title):
    if title == 'LFP':
        pot = get_extracellular(h, pop_names, time_pts, ele_pos, variable) #LOWPASS
        lfp = lowpassfilter(pot)
    else:
        lfp = get_extracellular(h, pop_names, time_pts, ele_pos, variable)
    ax,im = plot_potential(ax, lfp, 0.05, title)
    return ax,im

def set_axis(ax, letter=None):
    if letter is not None:
        ax.text(0.05, 1.025, letter, fontsize=12, weight='bold', transform=ax.transAxes)
    #ax.set_aspect('equal')
    return ax

def consolidated_figure(h_dict, pop_names, ele_pos, time_pts, fig):
    z_steps = 2 #4 for the plots 1 for colorbar
    height_ratios = [1 for i in range(z_steps)]
    #height_ratios.append(0.07) #height of colorbar
    gs = gridspec.GridSpec(z_steps, 3, height_ratios=height_ratios)
    #extracellular
    ax = plt.subplot(gs[0, 0])
    title = 'Extracellular'
    h = h_dict[title]
    ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, 'i',  title)
    set_axis(ax, 'A')
    ax.set_ylabel('Electrode number')
    #LFP
    ax = plt.subplot(gs[0, 1])
    title = 'LFP'
    h = h_dict[title]
    ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, 'i',  title)
    set_axis(ax, 'B')
    #Only Passive
    ax = plt.subplot(gs[0, 2])
    title = 'Passive only'
    h = h_dict[title]
    ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, 'i',  title)
    set_axis(ax, 'C')
    #Only Passive soma and axon
    ax = plt.subplot(gs[1, 0])
    title = 'Passive soma and axon'
    h = h_dict[title]
    ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, 'i',  title)
    set_axis(ax, 'D')
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Electrode number')
    #cax = plt.subplot(gs[2, 0])
    #cbar = plt.colorbar(im, cax=cax, orientation='horizontal', extend='both')
    #passive axon
    ax = plt.subplot(gs[1, 1])
    title = 'Passive axon'
    h = h_dict[title]
    ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, 'i',  title)
    set_axis(ax, 'E')
    ax.set_xlabel('Time (ms)')
    #cax2 = plt.subplot(gs[2, 1])
    #cbar2 = plt.colorbar(im, cax=cax2, orientation='horizontal', extend='both')
    #closed fast Na
    ax = plt.subplot(gs[1, 2])
    title = 'Closed fast Sodium'
    h = h_dict[title]
    ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, 'i',  title)
    set_axis(ax, 'F')
    ax.set_xlabel('Time (ms)')
    #cax3 = plt.subplot(gs[2, 2])
    #cbar3 = plt.colorbar(im, cax=cax3, orientation='horizontal', extend='both')
    return fig

num_cmpts = [74, 74, 59, 59, 59, 59, 61, 61, 50, 59, 59, 59]
cell_range = [0,1000,1050,1140,1230,1320,1560,2360,2560,3060,3160,3260,3360]
num_cells = np.diff(cell_range) / 10 #10% MODEL
total_cmpts = list(num_cmpts*num_cells)
pop_names = ['pyrRS23','pyrFRB23','bask23','axax23','LTS23', 
             'spinstel4', 'tuftIB5', 'tuftRS5', 'nontuftRS6', 
             'bask56', 'axax56', 'LTS56']
h24 = h5.File('/media/cchintaluri/PersonalBackup/data/hela_data/small_awake05/1/traub_syn.h5', 'r')
h25 = h5.File('/media/cchintaluri/PersonalBackup/data/hela_data/small_awake05/2/traub_syn.h5', 'r')
h26 = h5.File('/media/cchintaluri/PersonalBackup/data/hela_data/small_awake05/3/traub_syn.h5', 'r')
h27 = h5.File('/media/cchintaluri/PersonalBackup/data/hela_data/small_awake05/4/traub_syn.h5', 'r')
h28 = h5.File('/media/cchintaluri/PersonalBackup/data/hela_data/small_awake05/6/traub_syn.h5', 'r')
hs = [h24, h25, h26, h27, h28]
h_dict = {'Extracellular':h24,
          'LFP':h24,
          'Passive only':h25,
          'Passive soma and axon':h26,
          'Passive axon':h27,
          'Closed fast Sodium':h28}


max_val_dict = {'i':0.05, 'i_gaba_a': 0.05, 'i_cap':0.05, 
                'i_pas':0.05, 'i_k':0.5, 'i_na':0.5, 'i_ca':0.5, 
                'i_cat':0.05, 'i_ar':0.05, 'AMPA+NMDA':0.5}
title_dict = {'i':'LFP', 'i_gaba_a': 'GABA A', 'i_cap':'Capacitive', 
              'i_pas':'Passive', 'i_k':'Potassium', 'i_na':'Sodium', 
              'i_ca':'Calcium', 'i_cat':'CAT', 'i_ar':'AR', 'AMPA+NMDA':'AMPA+NMDA'}
num_ele = 28
ele_pos = place_electrodes_1D(num_ele)
time_pts = 6000
fig = plt.figure(figsize=(12,9))
fig = consolidated_figure(h_dict, pop_names, ele_pos, time_pts, fig)
[h.close() for h in hs] #close all files memory overflow
plt.tight_layout()
plt.savefig('fig4.png', dpi=300)
#plt.show()
