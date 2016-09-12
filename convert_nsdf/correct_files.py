import h5py as h5
import numpy as np
import pandas as pd
from time import strftime

def correct_amda_nmda(h, path): 
    #pass pop_name entry to be corrected
    print path
    loc = h[path]
    try:
        loc["i_NMDA"] = loc["i_ampa"]
        loc["i_AMPA"] = loc["i_nmda"]
        del loc["i_ampa"], loc["i_nmda"]
    except KeyError:
        print 'No NMDA/AMPA here', path
    return h

def delete_all_DS(h, pop_name):
    for dataset in h[pop_name].itervalues():
        for attr_key in dataset.attrs.iterkeys():
            dataset.attrs.__delitem__(attr_key)
    return h

def add_time_attributes(h, pop_name):
    #Add attributes field, dt, tunit, tstart, unit
    for dataset in h[pop_name].itervalues():
        field_name = str(dataset.name.lstrip(pop_name))
        dataset.attrs.create('field', field_name) #ADD FIELD attrb
        dataset.attrs.create('dt', 0.1)
        dataset.attrs.create('tunit', 'ms')
        dataset.attrs.create('tstart', 0.0)
        if field_name.find('i') != -1:
            dataset.attrs.create('unit', 'nA')
        else:
            dataset.attrs.create('unit', 'mV')
    return h

def tie_data_map(d_set, m_set, name, axis=0):
    d_set.dims[axis].label = name
    d_set.dims.create_scale(m_set, name)
    d_set.dims[axis].attach_scale(m_set)
    m_set.attrs.create('NAME', data='source')

if __name__ == '__main__':
    cell_range = [0,1000,1050,1140,1230,1320,1560,2360,2560,3060,3160,3260,3360,3460,3560]
    pop_names = ['pyrRS23','pyrFRB23','bask23','axax23','LTS23', 
                 'spinstel4', 'tuftIB5', 'tuftRS5', 'nontuftRS6', 
                 'bask56', 'axax56', 'LTS56', 'TCR', 'nRT']
    num_cmpts = [74, 74, 59, 59, 59, 59, 61, 61, 50, 59, 59, 59, 137, 59]
    #h = h5.File('testing.h5', 'a')

    h = h5.File('traub.h5', 'a')

    ##CORRECT NMDA NAME, DELETE DS
    prefix = '/data/uniform/'
    for pop in pop_names:
        pop_name = prefix+pop
        correct_amda_nmda(h, pop_name)
        delete_all_DS(h, pop_name)
        print 'DONE for :', pop_name

    ##ADD ADDITIONAL ATTRIBUTES!
    prefix_datasets = '/data/uniform/'
    for pop in pop_names:
        pop_name = prefix_datasets + pop    
        for dataset in h[pop_name].itervalues():
            field_name = str(dataset.name.rsplit('/',1)[1])
            dataset.attrs.create('field', field_name) #ADD FIELD attrb
            dataset.attrs.create('dt', 0.1)
            dataset.attrs.create('tunit', 'ms')
            dataset.attrs.create('tstart', 0.0)
            if field_name.find('i') != -1:
                dataset.attrs.create('unit', 'nA')
            else:
                dataset.attrs.create('unit', 'mV')
        print 'ADDED attrbs for dataset: ', dataset.name
             
    ##ADD ADDITIONAL ATTRIBUTES!
    prefix_datasets = '/data/event/'
    for pop in pop_names:
        pop_name = prefix_datasets + pop    
        for dataset in h[pop_name].itervalues():
            #field_name = str(dataset.name.rsplit('/',1)[1])
            dataset.attrs.create('field', 'spiketime') #ADD FIELD attrb
            dataset.attrs.create('unit', 'ms')
        print 'ADDED attrbs for dataset: ', dataset.name

    #Morphology!
    duh_ = 0 #counter to keep track of cells
    all_pops = range(14)
    for pop_idx in all_pops: #each population do
        num_cells = cell_range[pop_idx+1] - cell_range[pop_idx] 
        cmpts_per_cell = num_cmpts[pop_idx]
        unique_names = []
        for cell_idx in range(cell_range[pop_idx], cell_range[pop_idx+1]):
            for compt_n in range(cmpts_per_cell): #generating names for compartments
                unique_names.append(str(cell_idx)+'/'+str(compt_n+1))
        print 'Morphology for:', num_cells*cmpts_per_cell, pop_names[pop_idx]
        sp_type = np.dtype([('x0', np.float64),('y0', np.float64),('z0', np.float64),
                            ('x1', np.float64),('y1', np.float64),('z1', np.float64),
                            ('d', np.float64)])
        d_set = h.create_dataset('/data/static/morphology/'+ pop_names[pop_idx], 
                                shape=(num_cells*cmpts_per_cell,), dtype=sp_type)
        arr = pd.read_csv('punkty0.dat', sep='\t', usecols=[2,3,4,5,6,7,8], 
                          nrows=num_cells*cmpts_per_cell, skiprows= duh_,
                          names=['x0','y0','z0','x1','y1','z1','d'], delim_whitespace=True)
        d_set[0:num_cells*cmpts_per_cell, 'x0'] = arr.values[:,0] #something wrong with assigning directly
        d_set[0:num_cells*cmpts_per_cell, 'y0'] = arr.values[:,1]
        d_set[0:num_cells*cmpts_per_cell, 'z0'] = arr.values[:,2]
        d_set[0:num_cells*cmpts_per_cell, 'x1'] = arr.values[:,3]
        d_set[0:num_cells*cmpts_per_cell, 'y1'] = arr.values[:,4]
        d_set[0:num_cells*cmpts_per_cell, 'z1'] = arr.values[:,5]
        d_set[0:num_cells*cmpts_per_cell, 'd'] = arr.values[:,6]
        d_set.attrs.create('unit', ['um','um','um','um','um','um','um'])
        d_set.attrs.create('field', ['x0','y0','z0','x1','y1','z1','d'])
        m_set = h.create_dataset('/map/static/'+pop_names[pop_idx]+'_names', data=unique_names)
        tie_data_map(d_set, m_set, name='source', axis=0)
        duh_ += num_cells*cmpts_per_cell
    
    h.close()
