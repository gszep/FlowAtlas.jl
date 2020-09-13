from xml.etree.cElementTree import ElementTree
from networkx import DiGraph,draw,draw_networkx_edges,set_node_attributes
from pathlib import Path

from matplotlib.patches import Polygon,PathPatch
import fcsparser

from pandas import DataFrame,MultiIndex,SparseDtype,concat
from numpy import array,unique,roll,any,linspace
from re import search

import mpl_scatter_density
from matplotlib.scale import SymmetricalLogTransform
from matplotlib.colors import LinearSegmentedColormap,LogNorm

from matplotlib.pyplot import figure,show,subplot,axes,get_cmap,register_cmap
from networkx.drawing.nx_pydot import graphviz_layout

################################################# add transparency to colormaps
color_array = get_cmap('hot_r')(range(256))
color_array[:,-1] = linspace(1.0,0.0,256)[::-1]
register_cmap(cmap=LinearSegmentedColormap.from_list(name='hot_r',colors=color_array))

color_array = get_cmap('Blues')(range(256))
color_array[:,-1] = linspace(1.0,0.0,256)[::-1]
register_cmap(cmap=LinearSegmentedColormap.from_list(name='Blues',colors=color_array))

class Workspace(object):
	"""Parses FlowJo workspace file, metadata and gating strategy"""
	
	namespaces = {
		"xsi" : "http://www.w3.org/2001/XMLSchema-instance",
		"gating" : "http://www.isac-net.org/std/Gating-ML/v2.0/gating",
		"transformations" : "http://www.isac-net.org/std/Gating-ML/v2.0/transformations",
		"datatypes" : "http://www.isac-net.org/std/Gating-ML/v2.0/datatypes"
	}

	def __init__(self, path=str, channels={}, condition_parser = lambda x : {}):
		
		# parse samples from workspace
		parent_dir = str(Path(path).parent)
		samples,metadatas = [],[]

		# record dimensions for verification
		num_channels,num_labels = [],[]
		ncells = 0

		# format data for multicolumns
		for key,value in channels.items() :
			channels[key] = ('data',value)

		for sample in ElementTree(file=path).find('SampleList') :
	
			# truncate paths to parent folder level
			uri = sample.find('DataSet').get('uri')
			path_parts = list(Path(uri).parts)

			truncate_index = path_parts.index(parent_dir)
			uri = str(Path(*path_parts[truncate_index:]))

			############################################################ import data
			_,data = fcsparser.parse(path=uri,channel_naming='$PnN',dtype='float64')

			# get labels from gating
			gating = self.__gating__(sample)
			labels = self.apply( data, gating )

			####################################### rename channels according to provided map
			data.rename(columns=channels,inplace=True)
			set_node_attributes(gating,{ id: [ channels[key][-1] for key in gating.nodes[id]['channels'] ] for id in gating },'channels')

			metadata = DataFrame({'gating':[gating]})
			data = data.join(labels) # join labels and data to one frame
			data.columns = MultiIndex.from_tuples(data.columns)

			###################################### use conditions as multi-index
			for condition,value in condition_parser(uri).items() :
				data[condition] = value
				metadata[condition] = value

			data.index.name = 'cell'
			metadata.set_index( list(condition_parser(uri).keys()), inplace=True )
			data.set_index( list(condition_parser(uri).keys()), inplace=True, append=True )
			data.index = data.index.reorder_levels(roll(data.index.names,-1))

			metadatas.append(metadata)
			samples.append(data)

			######################################## verification stats
			cells,nchannels = data.shape
			cells,nlabels = labels.shape

			num_channels.append(nchannels)
			num_labels.append(nlabels)
			ncells += cells

		self.frame = concat(samples,copy=False)
		self.metadata = concat(metadatas,copy=False)

		############################## converting indexes to category types for efficiency
		for level in range(self.frame.index.nlevels-1) :
			self.frame.index.set_levels(
				self.frame.index.levels[level].astype('category'),
				level=level, inplace=True)

		for level in range(self.frame.columns.nlevels) :
			self.frame.columns.set_levels(
				self.frame.columns.levels[level].astype('category'),
				level=level, inplace=True)

		for level in range(self.metadata.index.nlevels) :
			self.metadata.index.set_levels(
				self.metadata.index.levels[level].astype('category'),
				level=level, inplace=True)

		self.frame.sort_index(inplace=True,level=range(self.frame.index.nlevels-1))
		self.metadata.sort_index(inplace=True)
		self.index = self.metadata.index

		############################################################################ verifications
		assert ncells == self.frame.shape[0], 'cells missing after concatination'
		assert len(unique(num_channels)) == 1, 'different number of channels per sample'

		assert len(unique(num_labels)) == 1, 'different number of populations per sample'
		assert self.frame.isnull().sum().sum() == 0, 'missing data after concatination; check channel names'


	def __gating__(self, sample) :
		'''parse gating strategy as networkx object'''

		# extract gating strategy
		populations = list(sample.iter('Population'))
		assert len(populations)>0, 'gating strategy not found for {}'.format(sample)
		graph = DiGraph()

		for population in populations :
			name = '__gate__'+population.get('name')
			Gate = population.find('Gate')

			id = Gate.get("{%s}id" % self.namespaces['gating'])
			parent_id = Gate.get("{%s}parent_id" % self.namespaces['gating'])

			for gate in Gate :
				assert 'PolygonGate' in str(gate), '{} not supported'.format(
					str(gate).replace('{'+self.namespaces['gating']+'}','gating:'))
				
				channels = []
				for dimension in gate.findall('gating:dimension/datatypes:fcs-dimension',self.namespaces) :
					channels.append( dimension.get("{%s}name" % self.namespaces['datatypes']) )

				vertexes = []
				for vertex in gate.findall('gating:vertex/gating:coordinate',self.namespaces) :
					vertexes.append( vertex.get("{%s}value" % self.namespaces['datatypes']) )

				vertexes = array(vertexes).astype(float).reshape(-1,2)
				graph.add_node( id, parent_id=parent_id, gate_name=name, channels=channels,
					gate = Polygon(vertexes,fill=False,edgecolor='midnightblue') )

				if parent_id is not None :
					graph.add_edge(parent_id,id)

		return graph


	def gate(self, data, gating, id):
		'''apply gate to DataFrame from node[id]'''

		node = gating.nodes[id]
		name,channels,gate = node['gate_name'],node['channels'],node['gate']

		# current gate 
		data[name] = gate.get_path().contains_points(data.filter(channels))

		if gating.in_degree(id)>0 : # and parents
			parent = array([ data[gating.nodes[parent]['gate_name']] for parent in gating.predecessors(id) ])
			data[name] &= any(parent,axis=0)

		# iterate through children
		for child in gating.successors(id) :
			data = self.gate(data,gating,child)

		return data


	def apply(self, data, gating):
		'''apply gating tree DataFrame and return boolean one-hot array'''

		roots = [ id for id in gating if gating.in_degree(id)==0 ]
		for id in roots :
			data = self.gate(data,gating,id)
			
		labels = data.filter([ gating.nodes[id]['gate_name'] for id in gating ])
		data.drop( columns = labels.columns, inplace=True)

		labels.columns = labels.columns.str.lstrip('__gate__')
		labels.columns = [ ('labels',column) for column in labels.columns ]
		return labels   #.astype(SparseDtype(bool, fill_value=False))

				
	def show(self, id, xlim=(-1e3,3e3), ylim=(-1e3,3e3), vmin=1, vmax=100,
			 linthresh=1e3, linscale=0.5, plot_size = 0.08):
		'''display gating heirarchy'''

		scale = SymmetricalLogTransform(base=10, linthresh=linthresh, linscale=linscale).transform
		norm = LogNorm(vmin=vmin,vmax=vmax)

		data,metadata =  self.frame.xs(id), self.metadata.xs(id)
		positions = graphviz_layout(metadata.gating, prog='dot')

		fig = figure(figsize=(30,30))
		fig.suptitle('file_index = {}'.format(id), fontsize=16)
		ax = subplot(111)

		transform =  fig.transFigure.inverted().transform
		draw_networkx_edges(metadata.gating, positions, ax=ax,
				edge_color='cornflowerblue', width=5, arrowsize=1)

		for id in metadata.gating:

			node = metadata.gating.nodes[id]
			xchannel,ychannel = node['channels']

			name = node['gate_name'].lstrip('__gate__')
			gate = node['gate']

			x,y = transform(ax.transData.transform(positions[id])) # axes coordinates
			gate_axes = axes([x-plot_size/2.0,y-plot_size/2.0,plot_size,plot_size],projection='scatter_density')

			gate_axes.set_aspect('equal')
			gate_axes.set_title(name)

			################################################################### gate
			gate.xy = scale(gate.xy).reshape(-1,2)
			gate_axes.add_patch(gate)

			################################################################### populations
			gate_mask = data.labels[name]
			if metadata.gating.in_degree(id)==0 : # whole population

				gate_axes.scatter(        scale(data.data[xchannel][~gate_mask]),scale(data.data[ychannel][~gate_mask]), color='gold', s=0.1, zorder=-1)
				if len(data.data[xchannel][~gate_mask]) > 0 :
					gate_axes.scatter_density(scale(data.data[xchannel][~gate_mask]),scale(data.data[ychannel][~gate_mask]), cmap='hot_r', norm=norm, zorder=1)

			else : # parent populations
				for parent in metadata.gating.predecessors(id) :
					mask = data.labels[metadata.gating.nodes[parent]['gate_name'].lstrip('__gate__')] & (~gate_mask)

					gate_axes.scatter(        scale(data.data[xchannel][mask]),scale(data.data[ychannel][mask]), color='gold', s=0.1, zorder=-1)
					if len(data.data[xchannel][mask]) > 0 :
						gate_axes.scatter_density(scale(data.data[xchannel][mask]),scale(data.data[ychannel][mask]), cmap='hot_r', norm=norm, zorder=1)

			# gated population
			gate_axes.scatter(        scale(data.data[xchannel][gate_mask]),scale(data.data[ychannel][gate_mask]), color='cornflowerblue', s=0.1, zorder=-1)
			if len(data.data[xchannel][gate_mask]) > 0 :
				gate_axes.scatter_density(scale(data.data[xchannel][gate_mask]),scale(data.data[ychannel][gate_mask]), cmap='Blues', norm=norm, zorder=1)

			gate_axes.set_xlim(*xlim)
			gate_axes.set_ylim(*ylim)

			gate_axes.set_xticks([])
			gate_axes.set_yticks([])

			gate_axes.set_xlabel(xchannel)
			gate_axes.set_ylabel(ychannel)

			gate_axes.xaxis.label.set_size(10)
			gate_axes.yaxis.label.set_size(10)

		ax.axis('off')
		show()