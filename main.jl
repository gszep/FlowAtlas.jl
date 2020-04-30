using GigaSOM
using Glob

data_path = "./data"
file_paths = glob( joinpath(data_path,"*","*filtered.fcs"))

params, fcsmatrix = loadFCS("data/423C/blood/Donor_Pre-Blood_023.fcs")  # load the FCS file
fcsmatrix

eval(:dataSet)
loadFCSSet(:dataSet,file_paths)

exprs = fcsmatrix[:,1:13]  # extract only the data columns with expression values

som = initGigaSOM(exprs, 20, 20)    # random initialization of the SOM codebook
som = trainGigaSOM(som, exprs)      # SOM training
clusters = mapToGigaSOM(som, exprs) # extraction of per-cell cluster IDs
e = embedGigaSOM(som, exprs)        # EmbedSOM projection to 2D
params

fcsmatrix
