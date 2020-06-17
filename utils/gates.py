from xml.etree.cElementTree import ElementTree
from networkx import DiGraph,draw,draw_networkx_edges,draw_networkx_labels,get_node_attributes,spring_layout

from FlowCytometryTools import PolyGate
from pandas import Series
from numpy import array

from matplotlib.pyplot import figure,show
from networkx.drawing.nx_pydot import graphviz_layout

class GatingSchema(object):
	"""Parses GatingML documents to PolyGates within networkx objects"""
	
	namespaces = {
		"xsi" : "http://www.w3.org/2001/XMLSchema-instance",
		"gating" : "http://www.isac-net.org/std/Gating-ML/v2.0/gating",
		"transformations" : "http://www.isac-net.org/std/Gating-ML/v2.0/transformations",
		"datatypes" : "http://www.isac-net.org/std/Gating-ML/v2.0/datatypes"
	}

	def __init__(self, path=str):
		
		self.graph = DiGraph()
		for gate in ElementTree(file=path).getroot() :

			assert 'PolygonGate' in str(gate), '{} not supported'.format(
				str(gate).replace('{'+self.namespaces['gating']+'}','gating:'))
			
			id = gate.get("{%s}id" % self.namespaces['gating'])
			gate_name = gate.get("{%s}name" % self.namespaces['gating'])
			parent_id = gate.get("{%s}parent_id" % self.namespaces['gating'])
			
			channels = []
			for dimension in gate.findall('gating:dimension/datatypes:fcs-dimension',self.namespaces) :
				channels.append( dimension.get("{%s}name" % self.namespaces['datatypes']) )

			vertexes = []
			for vertex in gate.findall('gating:vertex/gating:coordinate',self.namespaces) :
				vertexes.append( vertex.get("{%s}value" % self.namespaces['datatypes']) )

			vertexes = array(vertexes).astype(float).reshape(-1,2)
			vertexes = list(zip(*vertexes.T))

			self.graph.add_node( id, parent_id=parent_id,
				gate=PolyGate(vertexes,channels=channels, name=str(gate_name)), gate_name=gate_name )
			if parent_id is not None :
				self.graph.add_edge(parent_id,id)


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


	def get_labels(self, FCMeasurement):
		'''apply gating tree to FCMeasurement object and return row-level labels'''
		
		labels = Series( {'label':None}, index=FCMeasurement.data.index, dtype='string')
		population_ids = ( id for id in self.graph.nodes() if self.graph.out_degree(id)==0 and self.graph.in_degree(id)==1 )

		for id in population_ids :
			
			label = self.graph.nodes[id]['gate_name']
			label_index = self.apply(FCMeasurement,id).data.index
			labels.loc[label_index] = label

		return labels

				
	def show(self):
		'''display gating heirarchy'''

		positions = spring_layout(self.graph,k=0.2)

		figure(figsize=(6,6))
		draw(self.graph, positions, node_size=100,
				node_color='w',edgecolors='',linewidths=3)

		draw_networkx_edges(self.graph, positions, 
				edge_color='r', width=5, arrowsize=1)
		draw_networkx_labels(self.graph, positions, labels=get_node_attributes(self.graph,'gate_name'))
		show()