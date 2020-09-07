from warnings import filterwarnings
from matplotlib import MatplotlibDeprecationWarning
from anndata import ImplicitModificationWarning

filterwarnings('ignore', category=ImplicitModificationWarning)
filterwarnings('ignore', category=MatplotlibDeprecationWarning)

from progressbar import ProgressBar
from FlowCytometryTools import FCPlate
from .workspace import Workspace

from pandas import Series,DataFrame,Index
from scanpy import AnnData

from numpy.random import choice
from numpy import array

from glob import glob
from os.path import join
from re import search

def parser(file_path) :
    '''get tissue and patient id from filepath'''
    
    tissue = search('_(.+?)_',file_path).group()
    if 'Blood' in tissue : tissue = '_Blood_'

    patient = search('/[0-9]+C/',file_path).group()
    return (tissue[1:-1],patient[1:-1])


def load( path, patients, tissues, markers, channel_map={},
    ID='', pattern='*.fcs', subsample = None) :

    ################################################# load fcs files
    plate = FCPlate.from_dir(
        ID=ID, path=path, pattern=pattern,
        
        parser=parser, position_mapper= lambda x: x,
        row_labels=tissues,col_labels=patients)

    markers,types = zip(*markers)
    bar = ProgressBar(maxval=3*len(plate)).start()
    n = 0

    ################################################ load labels
    workspace_path = glob(join('data','*.wsp'))[0]
    workspace = Workspace(workspace_path,uri_parser=parser)

    labels = {}
    for idx in plate :
        labels[idx] = workspace.get_labels(plate[idx],idx)

        n += 1
        bar.update(n)

    ################################################ map channel names
    for idx in plate :
        
        plate[idx].data = plate[idx].data.rename(columns=channel_map)
        channel_differences = set(plate[idx].data.columns)-set(markers)

        if len(channel_differences) > 0 :
            print('Tissue: {}\tPatient: {}'.format(*idx))
            print(plate[idx].meta['_channels_'])
            raise ValueError('''{}\nChannel names not in markers. Change marker list or use channel_map'''.format(channel_differences))

        try :
            plate[idx].data = plate[idx].data.reindex(columns=markers)
            plate[idx].meta['_channel_names_'] = tuple(markers)

        except :
            print('Tissue: {}\tPatient: {}'.format(*idx))
            print(plate[idx].meta['_channels_'])
            raise ValueError('Channel renaming failed')

        if plate[idx].data.isna().any().any() :
            nan_channel = plate[idx].data.columns[plate[idx].data.isna().all().argmax()]
            raise ValueError('''{} contains NaNs in {}'''.format(idx,nan_channel))

        n += 1
        bar.update(n)

    ################################################# convert to AnnData object
    others = []
    for i,idx in enumerate(plate) :

        observations = DataFrame( index=plate[idx].data.index, columns=['tissue','patient','label'])
        observations['tissue'],observations['patient'] = idx
        observations['label'] = labels[idx]

        antigen = DataFrame(index=Index(plate[idx].data.columns,name='index'), columns=['type'])
        antigen['type'] = types

        annotated_data = AnnData( plate[idx].data, obs=observations, var=antigen)
        ################################################# subsample data
        if subsample != None :
            if subsample < annotated_data.n_obs :

                choices = choice( range(annotated_data.n_obs), size=subsample, replace=False )
                annotated_data = annotated_data[choices]

        if i == 1 : first = annotated_data
        else : others += [annotated_data]

        n += 1
        bar.update(n)

    dataset = first.concatenate(others)
    dataset.obs.drop(columns='batch',inplace=True)
    return dataset