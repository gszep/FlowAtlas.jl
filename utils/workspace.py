from xml.etree.cElementTree import ElementTree
from networkx import DiGraph,draw,draw_networkx_edges
from pathlib import Path

from FlowCytometryTools import PolyGate
from pandas import DataFrame
from numpy import array,arcsinh

from matplotlib.pyplot import figure,fill,show,title,subplot,axes
from networkx.drawing.nx_pydot import graphviz_layout

class Workspace(object):
	"""Parses FlowJo workspace file, metadata and gating strategy"""
	
	namespaces = {
		"xsi" : "http://www.w3.org/2001/XMLSchema-instance",
		"gating" : "http://www.isac-net.org/std/Gating-ML/v2.0/gating",
		"transformations" : "http://www.isac-net.org/std/Gating-ML/v2.0/transformations",
		"datatypes" : "http://www.isac-net.org/std/Gating-ML/v2.0/datatypes"
	}

	def __init__(self, path=str, uri_parser = None):
		
		# parse metadata from samples
		parent_dir = str(Path(path).parent)
		self.datasets = {}
		for sample in ElementTree(file=path).find('SampleList') :

			# parse channel and gating information
			channels = self.__channels__(sample)
			gating = self.__gating__(sample)
	
			# truncate paths to parent folder level
			uri = sample.find('DataSet').get('uri')
			path_parts = list(Path(uri).parts)

			truncate_index = path_parts.index(parent_dir)
			uri = str(Path(*path_parts[truncate_index:]))

			# record metadata per sample
			metadata = { 'uri': uri, 'channels': channels, 'gating': gating }

			if uri_parser == None : 
				self.datasets[uri] = metadata
			else : # custom parser
				self.datasets[uri_parser(uri)] = metadata

		self.ids = list(self.datasets.keys())


	def __channels__(self, sample) :
		'''parse channel-marker pairs with transformation lambdas'''
		parameters,channels = {},{}

		for keyword in sample.find('Keywords') :
			name = keyword.get('name')

			if name == '$PAR' : # get number of channels
				n_parameters = int(keyword.get('value'))

			elif '$P' == name[:2] : # extract channel/maker names
				parameters[name] = keyword.get('value')

		for n in range(1,n_parameters+1) : # match channel-marker pairs
			channels[parameters['$P'+str(n)+'N']] = { 'marker': parameters['$P'+str(n)+'S'] }

		# attach transfromation info
		return self.__transformations__(sample, channels)


	def __transformations__(self, sample, channels) :
		'''parse transformation as lambda functions'''

		for transformation in sample.find('Transformations') :
			channel = transformation.find("{%s}parameter"%self.namespaces['datatypes']).get('{%s}name'%self.namespaces['datatypes'])
			attributes = {}

			# parse transformation information
			name = transformation.tag.replace('{'+self.namespaces['transformations']+'}','')
			for key,value in transformation.attrib.items() :
				attributes[key.replace('{'+self.namespaces['transformations']+'}','')] = float(value)

			channels[channel]['transformation'] = self.transform(name,attributes)

		return channels


	def __gating__(self, sample) :
		'''parse gating strategy as networkx object'''

		# extract gating strategy
		populations = list(sample.iter('Population'))
		graph = DiGraph() if len(populations)>0 else None

		for population in populations :
			gate_name = population.get('name')
			Gate = population.find('Gate')

			id = Gate.get("{%s}id" % self.namespaces['gating'])
			parent_id = Gate.get("{%s}parent_id" % self.namespaces['gating'])

			for gate in Gate :
				assert 'PolygonGate' in str(gate), '{} not supported'.format(
					str(gate).replace('{'+self.namespaces['gating']+'}','gating:'))
				
				gate_channels = []
				for dimension in gate.findall('gating:dimension/datatypes:fcs-dimension',self.namespaces) :
					gate_channels.append( dimension.get("{%s}name" % self.namespaces['datatypes']) )

				vertexes = []
				for vertex in gate.findall('gating:vertex/gating:coordinate',self.namespaces) :
					vertexes.append( vertex.get("{%s}value" % self.namespaces['datatypes']) )
				
				vertexes = array(vertexes).astype(float).reshape(-1,2)
				vertexes = list(zip(*vertexes.T))

				graph.add_node( id, parent_id=parent_id,
					gate=PolyGate(vertexes, gate_channels, name=str(gate_name), region='in'), gate_name=gate_name )
				if parent_id is not None :
					graph.add_edge(parent_id,id)

		return graph


	def transform(self, name, attributes):
		'''return lamba function from attributes'''

		length = attributes['length']
		maxRange = attributes['maxRange']

		if name == 'fasinh' :
			assert maxRange == attributes['T']
			negative_decades = attributes['A']
			attributes['M'] # no idea what this parameter does

			width_basis = abs(attributes['W'])

		elif name == 'biex' :
			negative_decades = attributes['neg']
			attributes['pos'] # no idea what this parameter does
			width_basis = abs(attributes['width'])

		else :
			raise NotImplementedError('{} transformation not implemented'.format(name))
		
		return lambda x : arcsinh(array(x)/width_basis)


	def apply(self, FCMeasurement, id):
		'''apply gate to FCMeasurement object from node[id]'''

		node = self.graph.nodes[id]
		gate = node['gate']

		# match channel names between FCMeasurement object and GatingML file
		if '$PnS' in FCMeasurement.channels.columns and set(gate.channels).issubset(set(FCMeasurement.channels['$PnN'])):
			gate.channels = FCMeasurement.channels.set_index('$PnN').loc[gate.channels]['$PnS']

		# perform gate
		FCMeasurement = FCMeasurement.gate(gate)

		# continue performing gates if node has parents
		if self.graph.in_degree(id)==1 :
			return self.apply(FCMeasurement,node['parent_id'])
		else :
			return FCMeasurement


	def get_labels(self, FCMeasurement, sample_id):
		'''apply gating tree to FCMeasurement object and return row-level labels'''
		
		self.graph = self.datasets[sample_id]['gating']
		labels = DataFrame( {'label': [[]] }, index=FCMeasurement.data.index)

		if self.graph == None :
			return labels.label

		population_ids = ( id for id in self.graph if self.graph.out_degree(id)==0 and self.graph.in_degree(id)==1 )
		for id in population_ids :
			
			label = self.graph.nodes[id]['gate_name']
			label_index = self.apply(FCMeasurement.copy(),id).data.index
			labels.loc[label_index,'label'] = labels.loc[label_index,'label'].apply(lambda x: x+[label])
			
		return labels.label

				
	def show(self, sample_id, plot_size = 0.05):
		'''display gating heirarchy'''

		metadata = self.datasets[sample_id]
		positions = graphviz_layout(metadata['gating'], prog='dot')

		fig = figure(figsize=(20,20))
		fig.suptitle(metadata['uri'], fontsize=16)
		ax = subplot(111)

		transform =  fig.transFigure.inverted().transform
		draw_networkx_edges(metadata['gating'], positions, ax=ax,
				edge_color='r', width=5, arrowsize=1)

		for id in metadata['gating']:

			node = metadata['gating'].nodes[id]
			gate = node['gate']

			xchannel,ychannel = gate.channels
			xchannel = metadata['channels'][xchannel]
			ychannel = metadata['channels'][ychannel]

			x,y = transform(ax.transData.transform(positions[id])) # axes coordinates
			gate_axes = axes([x-plot_size/2.0,y-plot_size/2.0,plot_size,plot_size])
			gate_axes.set_aspect('equal')

			gate_axes.set_title(node['gate_name'])
			gate_axes.set_xlim(-4,8); gate_axes.set_ylim(-4,8)

			x,y = zip(*gate.vert) # display gate
			fill(xchannel['transformation'](x),ychannel['transformation'](y),fill=False,linewidth=3,color='orange')

			gate_axes.set_xlabel(xchannel['marker'])
			gate_axes.set_ylabel(ychannel['marker'])

			gate_axes.xaxis.label.set_size(10)
			gate_axes.yaxis.label.set_size(10)

			gate_axes.tick_params(axis='both', which='major', labelsize=6)

		ax.axis('off')
		show()