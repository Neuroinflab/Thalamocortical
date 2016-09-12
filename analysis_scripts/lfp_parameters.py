import numpy as np
import h5py as h5

cell_range = [0,1000,1050,1140,1230,1320,1560,2360,2560,3060,3160,3260,3360]
pop_names = ['pyrRS23','pyrFRB23','bask23','axax23','LTS23', 
             'spinstel4', 'tuftIB5', 'tuftRS5', 'nontuftRS6', 
             'bask56', 'axax56', 'LTS56']
num_cmpts = [74, 74, 59, 59, 59, 59, 61, 61, 50, 59, 59, 59]

def place_electrodes_1D(n):
    '''places n number of electrodes in ID along the column'''
    ele_x = np.ones((n,1))*25.
    ele_y = np.linspace(-2100., -80.0, num=n).reshape(n,1)
    ele_z = np.ones((n,1))*25.
    return np.hstack((ele_x, ele_y, ele_z))

def place_electrodes_2D(nx, ny):
    '''place nx*ny electrodes next to the column - like an MEA'''
    tot_ele = nx*ny
    zz = np.ones((tot_ele,1))*-25.
    xx, yy = np.mgrid[-375:375:np.complex(0,nx), -2250:450:np.complex(0,ny)]
    xx = xx.reshape(tot_ele, 1)
    yy = yy.reshape(tot_ele, 1)
    
    return np.hstack((xx, yy, zz))

def fetch_mid_pts(h, pop_name):
    all_pts = h['/data/static/morphology/'+pop_name]
    x = (all_pts['x0']+all_pts['x1']) / 2.
    y = (all_pts['y0']+all_pts['y1']) / 2.
    z = (all_pts['z0']+all_pts['z1']) / 2.
    x = x.reshape(x.size,1)
    y = y.reshape(y.size,1)
    z = z.reshape(z.size,1)
    return np.hstack((x, y, z))

# print place_electrodes_1D(10)
# print place_electrodes_2D(3, 3)
#h = h5.File('traub.h5', 'r')
h = h5.File('/home/cchintaluri/Hela_data_paper/repood/testing.h5', 'r')

ele_config = '1D'
ANIMATE = False

if ele_config == '1D':
    num_ele = 20
    ele_pos = place_electrodes_1D(num_ele)
    pot_filename = 'pot_sum_1D_'+str(num_ele)+'.npy'
elif ele_config == '2D':
    num_x, num_y = 16, 20
    num_ele = num_x*num_y
    ele_pos = place_electrodes_2D(num_x, num_y)
    pot_filename = 'pot_sum_2D_'+str(num_ele)+'.npy'
