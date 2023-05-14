from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import pkg_resources

#from ._0imports import *
import os
import sys

if sys.version_info[0] > 2:
    from urllib.request import urlretrieve
else:
    from urllib import urlretrieve
    input = raw_input

#import matplotlib
#if os.environ.get('DISPLAY', '') == '':
#    print('no display found. Using non-interactive Agg backend')
#    matplotlib.use('Agg')
#else:
#    matplotlib.use('TkAgg')

import glob
import time
import shutil
import pickle
import h5py
import astropy.units as u

__all__ = ["databases"]

class Database:

    def __init__(self, database_name, vital=False, date_to_update='daily', force_update=False, ask_size=None):

        self.database_name = database_name

        package_name = 'exotethys'
        info_file_name = '_0database.pickle'
        package_data = {package_name:'Passbands/*.pass'}
        last_update_file_name = 'database_last_update.txt'

        info_file_path = pkg_resources.resource_filename(package_name, info_file_name)
        package_path = os.path.join(os.path.expanduser('~'), '.{0}'.format(package_name))
        if not os.path.isdir(package_path):
            os.mkdir(package_path)

        self.package_path = package_path

        #self.directory_path = os.path.join(package_path, '{0}_{1}'.format(database_name, directory_name))
        self.directory_path = os.path.join(package_path, '{0}'.format(database_name))
        last_update_file_path = os.path.join(package_path, '{0}_{1}'.format(database_name, last_update_file_name))

        if date_to_update == 'daily':
            date_to_update = int(time.strftime('%y%m%d'))
        else:
            date_to_update = int(date_to_update)

        if os.path.isdir(self.directory_path):
            if force_update and len(glob.glob(os.path.join(self.directory_path, '*'))) != 0:
                shutil.rmtree(self.directory_path)
                os.mkdir(self.directory_path)	
                update = True
            elif len(glob.glob(os.path.join(self.directory_path, '*'))) == 0:
                update = True
            else:
                if not os.path.isfile(last_update_file_path):
                    update = True
                elif int(open(last_update_file_path).readlines()[0]) < date_to_update:
                    update = True
                else:
                    update = False
        else:
            os.mkdir(self.directory_path)
            update = True

        #with open(os.path.join(package_name, info_file_name), 'rb') as file:
        with open(info_file_path,'rb') as file: #Christophe correction
            dbx_files_dict = pickle.load(file)
            self.dbx_files = dbx_files_dict[database_name]

        #dbx_files = pickle.load(open(info_file_path, 'rb'))
        #dbx_files = dbx_files['{0}_{1}'.format(database_name, directory_name)]


#    def self_print(self):
#        print(self.database_name)
#        print(self.directory_name)
#        print(self.package_path)
#        print(self.directory_path)

    def get_file_content(self, dbx_file):
        abs_path_file = os.path.join(self.package_path, self.dbx_files[dbx_file]['local_path'])

        if not os.path.isfile(abs_path_file):
            print('Downloading... ', dbx_file)
            urlretrieve(self.dbx_files[dbx_file]['link'], os.path.join(self.package_path, self.dbx_files[dbx_file]['local_path']))
        else:
            print('File already here... ', dbx_file)
			
        file_ext = os.path.splitext(abs_path_file)[-1]
        if file_ext == '.pickle':
            try:
                with open(abs_path_file, 'rb') as file:
                    model_dict = pickle.load(file)
            except:
                print('Re-downloading... ', dbx_file)
                urlretrieve(self.dbx_files[dbx_file]['link'], os.path.join(self.package_path, self.dbx_files[dbx_file]['local_path']))
                with open(abs_path_file, 'rb') as file:
                    model_dict = pickle.load(file)
        elif file_ext in ['.h5', '.hdf5']:
            try:
                imodel_dict = h5py.File(abs_path_file,'r')
            except:
                print('Re-downloading... ', dbx_file)
                urlretrieve(self.dbx_files[dbx_file]['link'], os.path.join(self.package_path, self.dbx_files[dbx_file]['local_path']))
                imodel_dict = h5py.File(abs_path_file,'r')
            model_dict = {}
            for key in list(imodel_dict.keys()):
                model_dict[key] = imodel_dict[key][()] * u.Unit( imodel_dict[key].attrs['unit'] )
        return model_dict


	def get_filename_list(self):
		file_list = list(self.dbx_files.keys())
		return file_list


databases = {
"Phoenix_2012_13":Database('Phoenix_2012_13', date_to_update='200412', vital=True), "Phoenix_drift_2012":Database('Phoenix_drift_2012', date_to_update='201205', vital=True),
"Phoenix_2018":Database('Phoenix_2018', date_to_update='200412', vital=True), "Atlas_2000":Database('Atlas_2000', date_to_update='200403', vital=True), "Stagger_2015":Database('Stagger_2015', date_to_update='201001', vital=True), "Stagger_2018":Database('Stagger_2018', date_to_update='220216', vital=True), "MPS_Atlas_set1_2023":Database('MPS_Atlas_set1_2023', date_to_update='230512', vital=True)
}
