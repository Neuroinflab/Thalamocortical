##Meta data!!! 

dataset_locations = [
    '/512/200_awake0/traub.h5',
    '/512/100_awake0/traub.h5',
    '/512/50_awake0/traub.h5',
    '/512/25_awake0/traub.h5',
    '/512/12_5_awake0/traub.h5',
    '/512/8_awake0/traub.h5',
    '/512/4_awake0/traub.h5',
    '/512/2_awake0/traub.h5',
    '/512/kick/traub.h5',
    '/small_awake05/1/traub_syn.h5',
    '/small_awake05/2/traub_syn.h5',
    '/small_awake05/3/traub_syn.h5',
    '/small_awake05/4/traub_syn.h5',
    '/small_awake05/6/traub_syn.h5',
]

dataset_title = ['Oscillatory 200Hz no depol',
                 'Oscillatory 100Hz no depol', 
                 'Oscillatory 50Hz no depol', 
                 'Oscillatory 25Hz no depol', 
                 'Oscillatory 12.5Hz no depol', 
                 'Oscillatory 8Hz no depol',
                 'Oscillatory 4Hz no depol',
                 'Oscillatory 2Hz no depol',
                 'Whisker deflection',
                 'Whisker deflection 10% model',
                 'Whisker deflection passive currents',
                 'Whisker deflection passive soma & axon',
                 'Whisker deflection passive axon',
                 'Whisker deflection blocked Na current']


print len(dataset_title), len(dataset_locations)

import h5py as h5

def get_metadata(dataset_number):
    details_template = { 'title': dataset_title[dataset_number],
                         'description': dataset_description[dataset_number],
                         'creator' : 'Helena Glabska',
                         'contributor': 'Daniel Wojcik, Chaitanya Chintaluri',
                         'contact' : 'd.wojcik@nencki.gov.pl',
                         'created': strftime("%Y-%m-%d %H:%M:%S"),
                         'license' : 'ODbL 1.0',
                         'rights' : 'Copyright, Daniel K. Wojcik',
                         'tstart' : dataset_starttime[dataset_number],
                         'tend' : dataset_endtime[dataset_number],
                         'software' : 'Neuron 7.3, h5py',
                         'method' : 'CVODE, Crank-Nicholson, global dt',
                         'URL' : 'https://github.com/hglabska/Thalamocortical_imem',
                         'nsdf_version' : 1.0
                         }
    return details_template
                 

# with open("ODbLv1_0.txt", "r") as myfile:
#     license_ = myfile.read()
# str_type = h5.new_vlen(str)

all_datasets = range(14)

for dataset_number in all_datasets:
    path = '/media/cchintaluri/ed1076b7-6eda-45b5-a12b-ec0b1b69dafa/hela_data/'
    dataset_path = dataset_locations[dataset_number]
    h = h5.File(path + dataset_path, 'a')
    
    h.attrs.__delitem__('title')
    h.attrs.create('title', dataset_title[dataset_number])

    h.close()
    print 'Finished', dataset_path

#     #META DATA
#     metadata_dict = get_metadata(dataset_number)
#     for key, vals in metadata_dict.iteritems():
#         h.attrs.create(key, vals)

#     #LICENSE    
#     ds = h.create_dataset('/model/filecontents/LICENSE', shape=(1,), dtype=str_type)
#     ds[:] = license_

#     #URL
#     ds = h.create_dataset('/model/filecontents/url', shape=(1,), dtype=str_type)
#     ds[:] = 'https://github.com/hglabska/Thalamocortical_imem'

#     #README
#     ds = h.create_dataset('/model/filecontents/README', shape=(1,), dtype=str_type)
#     ds[:] = 'Copyright, Daniel K. Wojcik (d.wojcik@nencki.gov.pl)\n\n\n LICENSE: ODbL 1.0 license\n\n This dataset is made available under the Open Database License: http://opendatacommons.org/licenses/odbl/1.0/. Any rights in individual contents of the database are licensed under the Database Contents License: http://opendatacommons.org/licenses/dbcl/1.0/'

#     h.close()
#     print 'Finished', dataset_path
