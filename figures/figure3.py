import scipy.spatial
from scipy import signal
import numpy as np
import h5py as h5
import brewer2mpl
from matplotlib import rcParams, cm, gridspec
import matplotlib.pyplot as plt

rcParams.update({'font.size': 8})
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
    plt.title(title)
    #plt.xlim(xmin=2500, xmax=4500)
    #plt.colorbar(extend='both')
    return ax, im

def plot_raster(h, ax, title):
    set2 = brewer2mpl.get_map('Set3', 'qualitative', 12).mpl_colors
    set2.insert(0, 'crimson')
    set2.append('dodgerblue')
    cell_range = [0,1000,1050,1140,1230,1320,1560,2360,2560,3060,3160,3260,3360,3460,3560]
    pop_names = ['pyrRS23','pyrFRB23','bask23','axax23','LTS23', 
                 'spinstel4', 'tuftIB5', 'tuftRS5', 'nontuftRS6', 
                 'bask56', 'axax56', 'LTS56', 'TCR', 'nRT']
    inh_cells = ['bask23','axax23','LTS23', 
                 'bask56', 'axax56', 'LTS56', 'nRT']
    cells = np.array(cell_range)/ 10 #10% MODEL
    cells_spike = np.diff(cells)
    intra_cell_gap = 6
    intra_pop_gap = 15
    count = 0
    for ii,pop_name in enumerate(pop_names):
        cell_count = cells_spike[ii]
        all_spikes = h['/data/event/'+pop_name+'/spikes'].value
        int_spikes = all_spikes[0:cell_count]
        color = set2[ii]
        if pop_name in inh_cells:
            marker = 'v'
        else:
            marker = '^'
        for kk in range(cell_count):
            x = int_spikes[kk]
            y = np.zeros_like(x) + count
            #ax.plot(x, y, '.', color=color)
            ax.scatter(x, y, s=7, alpha=0.8,facecolor=color, 
                       marker=marker, linewidth='0.2', edgecolor='gray') #edgecolor='none'
            count -= intra_cell_gap
        count -= intra_pop_gap
    mid_pts = cells_spike / 2
    tick_pts = -1*(mid_pts + cells[:-1])*intra_cell_gap
    for ii in range(len(pop_names)):
        tick_pts[ii] -= (intra_pop_gap*ii) 
    ax.set_yticks(tick_pts)
    ax.set_yticklabels(pop_names)
    plt.plot(np.zeros((5000))+300, np.arange(-4500,500), color='k')
    plt.xlim((275, 425))
    plt.ylim((-2400 ,50))
    plt.xticks(np.arange(300, 500, 100), np.arange(300, 500, 100))
    plt.title(title)
    return ax

def get_specific(ax, h, pop_names, time_pts, ele_pos, variable, max_val_dict, title_dict):
    if variable == 'AMPA+NMDA':
        pot = get_extracellular(h, pop_names, time_pts, ele_pos, 'i_AMPA')
        pot += get_extracellular(h, pop_names, time_pts, ele_pos, 'i_NMDA')
    else:
        pot = get_extracellular(h, pop_names, time_pts, ele_pos, variable)
    lfp = lowpassfilter(pot)
    ax,im = plot_potential(ax, lfp, max_val_dict[variable], title_dict[variable])
    return ax,im

def consolidated_figure(h, pop_names, ele_pos, time_pts, fig, max_val_dict, title_dict):
    z_steps = 4 #4 for the plots 1 for colorbar
    height_ratios = [1 for i in range(z_steps)]
    height_ratios.append(0.07) #height of colorbar
    gs = gridspec.GridSpec(z_steps+1, 3, height_ratios=height_ratios)
    #Spikes
    ax = plt.subplot(gs[0, 0])
    spike_ax = plot_raster(h, ax, 'Spikes')
    #LFP
    ax = plt.subplot(gs[0, 1])
    variable = 'i'
    ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, variable, max_val_dict, title_dict)
    #AMPA + NMDA
    ax = plt.subplot(gs[0, 2])
    variable = 'AMPA+NMDA'
    ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, variable, max_val_dict, title_dict)
    #GABA A
    ax = plt.subplot(gs[1, 0])
    variable = 'i_gaba_a'
    ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, variable, max_val_dict, title_dict)
    plt.ylabel('Electrode depth ($\mu$m)')
    #I Cap
    ax = plt.subplot(gs[1, 1])
    variable = 'i_cap'
    ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, variable, max_val_dict, title_dict)
    #I K
    ax = plt.subplot(gs[1, 2])
    variable = 'i_k'
    ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, variable, max_val_dict, title_dict)
    #I Pas
    ax = plt.subplot(gs[2, 0])
    variable = 'i_pas'
    ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, variable, max_val_dict, title_dict)
    plt.ylabel('Electrode depth ($\mu$m)')
    #I Na
    ax = plt.subplot(gs[2, 2])
    variable = 'i_na'
    ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, variable, max_val_dict, title_dict)
    #I Cat
    ax = plt.subplot(gs[3, 0])
    variable = 'i_cat'
    ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, variable, max_val_dict, title_dict)
    plt.ylabel('Electrode depth ($\mu$m)')
    plt.xlabel('Time (ms)')
    cax1 = plt.subplot(gs[4, 0])
    cbar1 = plt.colorbar(im, cax=cax1, orientation='horizontal', extend='both')
    #I AR
    ax = plt.subplot(gs[3, 1])
    variable = 'i_ar'
    ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, variable, max_val_dict, title_dict)
    plt.xlabel('Time (ms)')
    cax2 = plt.subplot(gs[4, 1])
    cbar2 = plt.colorbar(im, cax=cax2, orientation='horizontal', extend='both')
    #I Ca
    ax = plt.subplot(gs[3, 2])
    variable = 'i_ca'
    ax, im = get_specific(ax, h, pop_names, time_pts, ele_pos, variable, max_val_dict, title_dict)
    plt.xlabel('Time (ms)')
    cax3 = plt.subplot(gs[4, 2])
    cbar3 = plt.colorbar(im, cax=cax3, orientation='horizontal', extend='both')
    return fig
    
num_cmpts = [74, 74, 59, 59, 59, 59, 61, 61, 50, 59, 59, 59]
cell_range = [0,1000,1050,1140,1230,1320,1560,2360,2560,3060,3160,3260,3360]
num_cells = np.diff(cell_range) / 10 #10% MODEL
total_cmpts = list(num_cmpts*num_cells)
pop_names = ['pyrRS23','pyrFRB23','bask23','axax23','LTS23', 
             'spinstel4', 'tuftIB5', 'tuftRS5', 'nontuftRS6', 
             'bask56', 'axax56', 'LTS56']
h = h5.File('/media/cchintaluri/PersonalBackup/data/hela_data/small_awake05/1/traub_syn.h5', 'r')

max_val_dict = {'i':0.05, 'i_gaba_a': 0.05, 'i_cap':0.05, 
                'i_pas':0.05, 'i_k':0.5, 'i_na':0.5, 'i_ca':0.5, 
                'i_cat':0.05, 'i_ar':0.05, 'AMPA+NMDA':0.5}
title_dict = {'i':'LFP', 'i_gaba_a': 'GABA A', 'i_cap':'Capacitive', 
              'i_pas':'Passive', 'i_k':'Potassium', 'i_na':'Sodium', 
              'i_ca':'Calcium', 'i_cat':'CAT', 'i_ar':'AR', 'AMPA+NMDA':'AMPA+NMDA'}
num_ele = 30
ele_pos = place_electrodes_1D(num_ele)
time_pts = 6000

#pot = get_extracellular(h, pop_names, time_pts, ele_pos, 'i')
#lfp = lowpassfilter(pot)
#fig = plt.figure()#figsize=(4,6))
#ax = subplot(111)
#ax1, im = plot_potential(ax, lfp, 0.05, 'LFP')
#ax1 = plot_raster(h, ax, 'Spikes')

fig = plt.figure(figsize=(12,18))
fig = consolidated_figure(h, pop_names, ele_pos, time_pts, fig, max_val_dict, title_dict)
plt.tight_layout()
plt.savefig('fig3.png', dpi=600)
#plt.show()
