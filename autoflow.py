from progressbar import ProgressBar
from utils.features import estimate_cutoff
from utils.fcs import load

from pandas import DataFrame
from numpy import arcsinh,mean


class AutoFlow:
    def __init__(self, ID, path, pattern, tissues, patients, markers,
                 channel_map = {}, subsample = 3000, cofactor = 250.0):

        ###################### load using FlowCytometryTools
        self.dataset = load( ID = ID, path = path, pattern = pattern,
            tissues = tissues, patients = patients, markers = markers,
            channel_map = channel_map, subsample=subsample
        )

        ###################### apply biexponential transfromation
        self.dataset.X = arcsinh(self.dataset.X/cofactor)

        ###################### store unprocessed dataset
        self.unprocessed = self.dataset.obs.join(
            DataFrame(self.dataset.X, index=self.dataset.obs.index,
                    columns=self.dataset.var.index))

    def batch_normalise(self):
        '''method uses controls and gaussian mixture models to detect
        positive and negative fluorescence regions.'''

        ################################################# batch normalisation
        batches = self.dataset.obs.reset_index().groupby(['patient','tissue'])
        bar = ProgressBar(maxval=len(batches)).start()
        n = 0

        for _,batch in batches :
            for j,channel in enumerate(self.dataset.X[batch.index].T) :
                
                # use gaussian mixture models to estimate +/- feature bounadary
                feature = channel - estimate_cutoff(channel)
                
                # scale negative and positive features independently to +/- 1
                feature[feature>0] /=  mean(feature[feature>0])
                feature[feature<0] /= -mean(feature[feature<0])
                
                self.dataset.X[batch.index,j] = feature
                
            n += 1
            bar.update(n)