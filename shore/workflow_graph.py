import plotly.graph_objects as go
import networkx as nx

class workflow:
    def __init__(self):
        self.graph = nx.Graph()
        self.id = 0

    def add_node(self, node, layer=None, state='default'):
        # Assign layer and state as attributes to the node
        self.graph.add_node(node, layer=layer, data='test', state=state)

    def add_edge(self, node1, node2):
        self.graph.add_edge(node1, node2)

    def add_instance(self, node1=None, node2=None, layer=None, state='default'):
        self.id += 1
        self.graph.add_node(node2, data=f'test{self.id}', layer=layer, state=state)
        self.add_edge(node1, node2)

    def show(self):
        """Visualize the multi-layered graph using Plotly."""
        pos = self.vertical_layered_layout()

        # Create edge traces
        edge_x = []
        edge_y = []

        for edge in self.graph.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.append(x0)
            edge_x.append(x1)
            edge_x.append(None)  # None separates edges
            edge_y.append(y0)
            edge_y.append(y1)
            edge_y.append(None)  # None separates edges

        # Create node traces
        node_x = []
        node_y = []
        node_colors = []

        # Define a color mapping based on states
        color_map = {
            'active': 'green',
            'inactive': 'red',
            'pending': 'yellow',
            'default': '#6175c1'  # Default color
        }

        node_text=[]
        for node in self.graph.nodes():
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)

            # Get the state of the node and assign a color
            state = self.graph.nodes[node]['state']
            color = color_map.get(state, '#6175c1')  # Fallback to default color if state not found
            node_colors.append(color)
            if self.graph.nodes[node]['layer'] in ['input','structure','xas results']:
                info_text = node
            elif self.graph.nodes[node]['layer']=='xas results':
                info_text=node.split('-')[1:]
            else:
                info_text = self.graph.nodes[node]['layer']
            node_text.append(info_text) 

        annotations = []
        for edge in self.graph.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]

            annotations.append(dict(
                ax=x0, ay=y0,
                axref='x', ayref='y',
                x=x1, y=y1,
                xref='x', yref='y',
                showarrow=True,
                arrowhead=2,
                arrowsize=5,
                arrowcolor='gray'
            ))

        # Create figure
        fig = go.Figure()

        # Add edges to figure
        fig.add_trace(go.Scatter(x=edge_x, y=edge_y,
                                 line=dict(width=0.5, color='#888'),
                                 mode='lines'))

        # Add nodes to figure with colors based on their states
        fig.add_trace(go.Scatter(x=node_x, y=node_y,
                                 mode='markers+text',
                                 text=node_text,
                                 textposition="top center",
                                 hoverinfo='text',
                                 
                                 marker=dict(
                                     size=30,
                                     color=node_colors,
                                     line_width=2),
                                
        ))
        
        xb,yb,xe,ye=self.cluster_layout()
            
        fig.add_shape(type="rect",
                              x0=xb - 1,
                              y0=yb + 1,
                              x1=xe + 1,
                              y1=ye - 1,
                              line=dict(color="RoyalBlue", width=0.1),
                              fillcolor="LightSkyBlue",
                              layer="below",
                              )
        # fig.add_annotation(
        #             x=xb,
        #             y=yb+0.2,
        #             text='cluster',
        #             showarrow=False,
        #             font=dict(color="black"),
        #             align="center"
        #         )


        fig.update_layout(showlegend=False,
                          hovermode='closest',
                          annotations=annotations,
                          margin=dict(b=0, l=0, r=0, t=0),
                          xaxis=dict(showgrid=False,
                                     zeroline=False,
                                     showticklabels=False),
                          yaxis=dict(showgrid=False,
                                     zeroline=False,
                                     showticklabels=False))

        fig.show()

    def vertical_layered_layout(self):
        """Create a layout for a vertically aligned multi-layered graph."""
        pos = {}
        layers = {}

        # Group nodes by layer
        for node, data in self.graph.nodes(data=True):
            layer = data['layer']
            if layer not in layers:
                layers[layer] = []
            layers[layer].append(node)

        # Assign positions in a vertical layout
        for layer_index, layer_nodes in enumerate(layers.values()):
            x = layer_index * 3  # Adjust spacing between layers
            for i, node in enumerate(layer_nodes):
                y = i * 3  # Adjust spacing between nodes in the same layer
                pos[node] = (x, y)

        return pos
    
    def cluster_layout(self):
        """Create a layout for a vertically aligned multi-layered graph."""
        pos = {}
        layers = {}

        # Group nodes by layer
        for node, data in self.graph.nodes(data=True):
            layer = data['layer']
            if layer not in layers:
                layers[layer] = []
            layers[layer].append(node)

        # Assign positions in a vertical layout
        for layer_index,(layer_name, layer_nodes) in enumerate(layers.items()):
            x = layer_index * 3  # Adjust spacing between layers
            if 'parsing' in layer_name:
                ytmp=[]
                for i, node in enumerate(layer_nodes):
                    y = i * 3  # Adjust spacing between nodes in the same layer
                    ytmp.append(y)
                xb=x
                yb=max(ytmp)
            if 'bse' in layer_name:
                ytmp=[]
                for i, node in enumerate(layer_nodes):
                    y = i * 3  # Adjust spacing between nodes in the same layer
                    ytmp.append(y)
                xe=x
                ye=min(ytmp)

        return xb,yb,xe,ye
