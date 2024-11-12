"""
A.Geondzhian, N. Davydov
"""
import logging
import numpy as np

import plotly.graph_objects as go
import plotly.express as px
import seaborn as sns
# from siman.core.structure import Structure


class visual:
    """
    class for crystal structure visualization
    """

    def __init__(self, atoms):
        """
        Initialize the visual
        Args:
        atoms  - for now only siman object which defines structure.
        """

        self.atoms = atoms
        self.param = dict(scatter=True,
                          projection=False,
                          projection_type='xy',
                          cell_vectors=True,
                          bonds=True,
                          bonds_color=False,
                          bonds_color_threshold=0.06,
                          bonds_elements=['Ni', 'O'],
                          bonds_length=True,
                          bonds_color_scale='Plasma',
                          bonds_threshold=2.4,
                          resolution=20,
                          unitcell=None,
                          )

        # self.fig = go.Figure()
        self.fig = None
        self.annotations = []
        self.bond_pairs = {}
        self.data = []
        self.abc = {0: "a", 1: "b", 2: 'c'}
        # test
        self.color_labels = {
            'Ni': {'color': 'rgb(77,0,242)', 'size': 0.3, 'opacity': 0.95},
            'O': {'color': 'rgb(235,0,0)', 'size': 0.15, 'opacity': 0.95},
            'Li': {'color': 'grey', 'size': 0.2, 'opacity': 0.5}
        }
        self.list_of_atoms=['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb ', 'Bi ', 'Po ', 'At ', 'Rn ', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U']
        self.list_of_colors=self.__get_color_palette('deep', len(self.list_of_atoms))
        for i,(item,icolor) in enumerate(zip(self.list_of_atoms,self.list_of_colors)):
            
            self.color_labels.update({
                item:dict(color=icolor, size=i*0.01, opacity=0.8)
            })
        self.color_labels.update( {
            'Ni': {'color': 'rgb(77,0,242)', 'size': 0.3, 'opacity': 0.95},
            'O': {'color': 'rgb(235,0,0)', 'size': 0.15, 'opacity': 0.95},
            'Li': {'color': 'grey', 'size': 0.2, 'opacity': 0.5}
        })


    def __get_color_palette(self,color_palette, num_colors):
        """Generate a list of RGB colors from a seaborn color palette."""
        palette = sns.color_palette(color_palette, num_colors)
        rgb_colors = ['rgb({},{},{})'.format(int(r * 255), int(g * 255), int(b * 255)) for r, g, b in palette]
        return rgb_colors


    def plot(self, param=None):
        """Method to create and fill with content plotly figure. 

        Args:
            param (dict, optional): dictionary of optional parameters to adjust visualizatoin. 
            Defaults to None.
        """
        if param is not None:
            self.param.update(param)

        for index, (name, pos) in enumerate(zip(self.atoms.get_chemical_symbols(), self.atoms.positions)):
            self.data.append(self.create_atom_instance(index, name, pos))
        # print(self.data.append)data=self.data
        self.fig = go.Figure(data=self.data)

        self._add_bonds()
        self._add_cell()
        self._general_fomarting()
        self.fig.data = self.fig.data[::-1]

    def plot2d(self, index, name, pos):
        """method to draw atoms as circiles in 2d

        Args:
            index (int): index which will be reflected in the hoverinfo
            name (str): atomic symbol
            pos (list): 2d position of the atom

        Returns:
            plotly graphical object: instance for the given atom
        """
        return go.Scatter(x=[pos[0]],
                          y=[pos[1]],
                          opacity=self.color_labels[name]['opacity'],
                          name=f'{name} index {index}',
                          marker=dict(size=self.color_labels[name]['size']*100,
                                      color=self.color_labels[name]['color']),
                          showlegend=False)

    def plot3d(self, index, name, pos):
        """method to draw atoms spheres with _ms or points with go.Scatter3d

        Args:
            index (int): index which will be reflected in the hoverinfo
            name (str): atomic symbol
            pos (list): 2d position of the atom

        Returns:
            plotly graphical object: instance for the given atom
        """
        if self.param['scatter']:
            return go.Scatter3d(x=[pos[0]],
                                y=[pos[1]],
                                z=[pos[2]],
                                opacity=self.color_labels[name]['opacity'],
                                name=f'{name} index {index}',
                                marker=dict(size=self.color_labels[name]['size']*50,
                                            color=self.color_labels[name]['color']),
                                showlegend=False)
        else:
            (x_pns_surface, y_pns_surface, z_pns_suraface) = self._ms(
                pos, self.color_labels[name]['size'])
            return go.Surface(x=x_pns_surface, y=y_pns_surface, z=z_pns_suraface,
                              opacity=self.color_labels[name]['opacity'],
                              name=f'{name} index {index}',
                              colorscale=[[0, self.color_labels[name]['color']],
                                          [1, self.color_labels[name]['color']]],
                              showscale=False)

    def apply_projection(self, pos: list):
        """ Method to process projections.

        Args:
            index (int): index of a given atom
            name (str): atomic symbol
            pos (list): list of x,y,z coordinates in agnstrom 
        """
        if self.param['projection']:
            if 'xz' in self.param['projection_type'].lower():  # type: ignore
                projected_x, projected_y = pos[0], pos[2]
            elif 'yz' in self.param['projection_type'].lower():  # type: ignore
                projected_x, projected_y = pos[1], pos[2]
            elif 'xy' in self.param['projection_type'].lower():  # type: ignore
                projected_x, projected_y = pos[0], pos[1]
            else:
                logging.warning(
                    'Something went wrong in apply_projection, using "xy" instead.')
                projected_x, projected_y = pos[0], pos[1]
            return [projected_x, projected_y]
        else:
            return pos

    def create_atom_instance(self, index, name, pos):
        """wrapper on plot3d plot2d

        Args:
            index (int): index of a given atom
            name (str): atomic symbol
            pos (list): postion in angst for 2d or 3d cases
        """
        if not self.param['projection']:
            return self.plot3d(index, name, pos)
        else:
            pos2d = self.apply_projection(pos)
            return self.plot2d(index, name, pos2d)

    def get_bonds(self):
        """Method to get information about bond. Requires some adjustment.
        """
        for i, (name, position) in enumerate(zip(self.atoms.get_chemical_symbols(), self.atoms.positions)):
            if name in self.param['bonds_elements']:  # type: ignore
                for j, (second_name, second_position) in enumerate(zip(self.atoms.get_chemical_symbols()[i:], self.atoms.positions[i:])):
                    # type: ignore
                    if second_name in self.param['bonds_elements']:
                        bond_length = np.sqrt(
                            np.dot(position-second_position, position-second_position))
                        if bond_length < self.param['bonds_threshold'] and bond_length > 0.:
                            self.bond_pairs.update({f'{i}{j}':
                                                    {'atom1': name,
                                                     'atom2': second_name,
                                                     'bond_length': bond_length,
                                                     'begin':  position,
                                                     'end': second_position,
                                                     'color': 'green',
                                                     'annotation':
                                                     {'x': ((position+second_position)/2)[0],
                                                      'y': ((position+second_position)/2)[1],
                                                         'z': ((position+second_position)/2)[2],
                                                      }
                                                     }})

    def bond_color(self,):
        """ Method to get color scale and mark bonds with diffrenet lengths. 
        """
        scale = [self.bond_pairs[k]['bond_length'] for k in self.bond_pairs]
        self.bond_pairs = {k: v for k, v in sorted(
            self.bond_pairs.items(), key=lambda item: item[1]['bond_length'])}
        edge = min(scale)
        tmp = int((-min(scale)+max(scale))/self.param['bonds_color_threshold'])
        n_colors = tmp if tmp != 0 else 1
        if n_colors != 1:
            colors = px.colors.sample_colorscale(self.param['bonds_color_scale'], [
                                                 n/(n_colors - 1) for n in range(n_colors)])
            step = (-min(scale)+max(scale))/(n_colors-1)
        else:
            colors = ['green']
            step = max(scale)

        index = 0
        for k in self.bond_pairs.keys():
            # print(index)
            if self.bond_pairs[k]['bond_length'] < edge+step:
                self.bond_pairs[k]['color'] = colors[index]
            else:
                edge = edge+step
                index += 1
                self.bond_pairs[k]['color'] = colors[index]

    def _general_fomarting(self):
        """Options for general formatting of plotly figure.
        """
        self.fig.update_layout(coloraxis_showscale=False)

        if not self.param['projection']:
            self.fig.update_scenes(xaxis_visible=False,  # type: ignore
                                   yaxis_visible=False, zaxis_visible=False)
            self.fig.update_layout(scene=dict(bgcolor='white'),
                                   plot_bgcolor='rgb(0,0,0)',
                                   paper_bgcolor='rgb(0,0,0)',
                                   margin=dict(r=0, b=0, l=0, t=0),
                                   width=1000,)
            self.fig.update_layout(scene_camera=dict(up=dict(x=0, y=1., z=0),
                                                     eye=dict(x=0, y=0., z=10)),
                                   margin=dict(r=0, b=0, l=0, t=0),)
            self.fig.update_scenes(camera_projection_type='orthographic')
            self.fig.update_layout(scene=dict(
                annotations=self.annotations, aspectmode='data'))
        else:
            # self.fig.update_scenes(xaxis_visible=False,  # type: ignore
            #                        yaxis_visible=False, yaxis_showticklabels=False, xaxis_showticklabels=False)
            self.fig.update_layout(plot_bgcolor='rgb(0,0,0,0)',
                                   paper_bgcolor='rgb(0,0,0,0)',
                                   margin=dict(r=50, b=50, l=50, t=50),
                                   hovermode='closest',
                                   )
            self.fig.update_layout(yaxis_scaleanchor="x",)
            self.fig.update_xaxes(visible=False)
            self.fig.update_yaxes(visible=False)
            self.fig.update_layout(showlegend=False)
            self.fig.update_layout(legend = dict(bgcolor = 'white'))
            # self.fig.update_layout(yaxis_visible=False, yaxis_showticklabels=False)
        self.fig.update_layout(coloraxis={'showscale': False})

    def _plot_bonds(self):
        """Internal method to draw bonds for 2d and 3d configurations.
        """
        for bond in self.bond_pairs:
            if not self.param['projection']:
                self.fig.add_trace(self._bond_line_3d(self.bond_pairs[bond]))
                if self.param['bonds_length']:
                    self.annotations.append(dict(
                        x=self.bond_pairs[bond]['annotation']['x'],
                        y=self.bond_pairs[bond]['annotation']['y'],
                        z=self.bond_pairs[bond]['annotation']['z'],
                        text=f"{np.round(self.bond_pairs[bond]['bond_length'],2)}",
                        showarrow=False))
            else:
                self.fig.add_trace(self._bond_line_2d(self.bond_pairs[bond]))
                if self.param['bonds_length']:
                    projected = self.apply_projection([self.bond_pairs[bond]['annotation']['x'],
                                                       self.bond_pairs[bond]['annotation']['y'],
                                                       self.bond_pairs[bond]['annotation']['z']])
                    self.fig.add_annotation(
                        x=projected[0],
                        y=projected[1],
                        text=f"{np.round(self.bond_pairs[bond]['bond_length'],2)}",
                        showarrow=False)

    def _add_bonds(self):
        """Internal method to add bonds on the figure. 
        """
        if self.param['bonds']:
            self.get_bonds()
            if self.param['bonds_color']:
                self.bond_color()
            self._plot_bonds()

    # def _add_cell(self):
    #     """Internal method to add cell vecotors on the figure.
    #     """
    #     if self.param['cell_vectors']:
    #         for index, vector in enumerate(self.atoms.cell):
    #             if not self.param['projection']:
    #                 self.fig.add_trace(self._cell_vectors_3d(vector, index))
    #             else:
    #                 self.fig.add_trace(self._cell_vectors_2d(vector, index))
    def _add_cell(self):
        """Internal method to add cell vectors on the figure as a 3D unit cell.
        """
        if self.param['cell_vectors']:
            # Get the cell vectors
            if self.param['unitcell']:
                cell_vectors = self.param['unitcell'].cell
            else:
                cell_vectors = self.atoms.cell
            unit=1
            # Define unit vectors for X, Y, Z directions
            x_vector = np.array([1, 0, 0])*unit  # Unit vector in X direction
            y_vector = np.array([0, 1, 0])*unit # Unit vector in Y direction
            z_vector = np.array([0, 0, 1])*unit # Unit vector in Z direction
            arrow_length = 1.
            # Add these vectors to the vertices for visualization
            origin = np.array([-2, -2, -2])
            
            # Calculate vertices of the unit cell
            vertices = np.array([
                [0, 0, 0],  # Origin
                cell_vectors[0],  # a
                cell_vectors[1],  # b
                cell_vectors[2],  # c
                cell_vectors[0] + cell_vectors[1],  # a + b
                cell_vectors[0] + cell_vectors[2],  # a + c
                cell_vectors[1] + cell_vectors[2],  # b + c
                cell_vectors.sum(axis=0)  # a + b + c
            ])

            if self.param['projection']:
                # Create a 3D unit cell representation
                self.fig.add_trace(self._create_unit_cell_2d(vertices))
            else:
                self.fig.add_trace(self._draw_boundary_lines(vertices))
                # If projection is enabled, you might want to handle that separately.
                self.fig.add_trace(go.Scatter3d(
                    x=[origin[0], origin[0] + x_vector[0]], 
                    y=[origin[1], origin[1] + x_vector[1]], 
                    z=[origin[2], origin[2] + x_vector[2]],
                    mode='lines',
                    line=dict(color='blue', width=5),
                    name='X Direction'
                ))

                self.fig.add_trace(go.Scatter3d(
                    x=[origin[0], origin[0] + y_vector[0]], 
                    y=[origin[1], origin[1] + y_vector[1]], 
                    z=[origin[2], origin[2] + y_vector[2]],
                    mode='lines',
                    line=dict(color='green', width=5),
                    name='Y Direction'
                ))

                self.fig.add_trace(go.Scatter3d(
                    x=[origin[0], origin[0] + z_vector[0]], 
                    y=[origin[1], origin[1] + z_vector[1]], 
                    z=[origin[2], origin[2] + z_vector[2]],
                    mode='lines',
                    line=dict(color='orange', width=5),
                    name='Z Direction'
                ))
                # X Direction Arrow
                self.fig.add_trace(go.Cone(
                    x=[origin[0]+1.], 
                    y=[origin[1]], 
                    z=[origin[2]], 
                    u=[arrow_length * x_vector[0]], 
                    v=[arrow_length * x_vector[1]], 
                    w=[arrow_length * x_vector[2]],
                    sizemode="absolute",
                    # size=0.1,
                    anchor="tail",
                    colorscale=[[0, 'blue'], [1, 'blue']],
                    showscale=False,
                    name='X Direction'
                ))

                # Y Direction Arrow
                self.fig.add_trace(go.Cone(
                    x=[origin[0]], 
                    y=[origin[1]+1.], 
                    z=[origin[2]], 
                    u=[arrow_length * y_vector[0]], 
                    v=[arrow_length * y_vector[1]], 
                    w=[arrow_length * y_vector[2]],
                    sizemode="absolute",
                    # size=0.1,
                    anchor="tail",
                    colorscale=[[0, 'green'], [1, 'green']],
                    showscale=False,
                    name='Y Direction'
                ))

                # Z Direction Arrow
                self.fig.add_trace(go.Cone(
                    x=[origin[0]], 
                    y=[origin[1]], 
                    z=[origin[2]+1.], 
                    u=[arrow_length * z_vector[0]], 
                    v=[arrow_length * z_vector[1]], 
                    w=[arrow_length * z_vector[2]],
                    sizemode="absolute",
                    # size=0.1,
                    anchor="tail",
                    colorscale=[[0, 'orange'], [1, 'orange']],
                    showscale=False,
                    name='Z Direction'
                ))
                

    def _draw_boundary_lines(self, vertices):
        """Draw boundary lines for the unit cell in 3D.
        
        Parameters:
        - vertices (ndarray): The vertices of the unit cell.
        
        Returns:
        - trace: A Plotly scatter3d trace representing the boundary lines.
        """
        # Define line segments connecting vertices to form edges of the parallelepiped
        edges = [
            [vertices[0], vertices[1]],  # Origin to a
            [vertices[0], vertices[2]],  # Origin to b
            [vertices[0], vertices[3]],  # Origin to c
            [vertices[1], vertices[4]],  # a to a + b
            [vertices[1], vertices[5]],  # a to a + c
            [vertices[2], vertices[4]],  # b to a + b
            [vertices[2], vertices[6]],  # b to b + c
            [vertices[3], vertices[5]],  # c to a + c
            [vertices[3], vertices[6]],  # c to b + c
            [vertices[4], vertices[7]],  # a + b to a + b + c
            [vertices[5], vertices[7]],  # a + c to a + b + c
            [vertices[6], vertices[7]]   # b + c to a + b + c
        ]

        x_lines = []
        y_lines = []
        z_lines = []

        for edge in edges:
            x_lines.extend([edge[0][0], edge[1][0], None])  # None is used to break lines in Plotly
            y_lines.extend([edge[0][1], edge[1][1], None])
            z_lines.extend([edge[0][2], edge[1][2], None])

        return dict(
            type='scatter3d',
            x=x_lines,
            y=y_lines,
            z=z_lines,
            mode='lines',
            line=dict(color='black', width=2),  showlegend=False, # Adjust color and width as needed
        )
    def _create_unit_cell_3d(self, vertices):
        """Create a mesh for the unit cell in 3D.
        
        Parameters:
        - vertices (ndarray): The vertices of the unit cell.
        
        Returns:
        - trace: A Plotly mesh3d trace representing the unit cell.
        """
        i = [0, 1, 1, 0, 0, 4, 4, 5, 5, 6, 6, 2, 2, 3]
        j = [1, 2, 4, 4, 5, 5, 6, 6, 3, 3, 1, 2]
        k = [2, 0, 5, 6, 7, 7, 4]

        return go.Mesh3d(
            x=vertices[:,0],
            y=vertices[:,1],
            z=vertices[:,2],
            i=i,
            j=j,
            k=k,
            opacity=0.5,
            color='rgba(100,100,250,0.5)', # Adjust color and opacity as needed
           
        )

    def _create_unit_cell_2d(self, vertices):
        """Create a projection of the unit cell in 2D.
        
        Parameters:
        - vertices (ndarray): The vertices of the unit cell.
        
        Returns:
        - trace: A Plotly scatter trace representing the projection of the unit cell.
        """
        # Project onto XY plane (assuming Z is up)
        projected_vertices = vertices[:, :2]  

        return dict(
            type='scatter',
            x=projected_vertices[:,0],
            y=projected_vertices[:,1],
            mode='lines',
            line=dict(color='blue', width=2),
            fill='toself',
            fillcolor='rgba(100,100,250,0.5)', # Adjust color and opacity as needed
        )

    def _cell_vectors_3d(self, vector, index):
        """Internal method to handle 3d unitcell vectors.

        Args:
            vector (list): end point for radius vector
            index (int): 0-a,1-b,2-c

        Returns:
            plotly.go: line and name of the vector
        """
        self.annotations.append(dict(x=vector[0]+0.5,
                                     y=vector[1]+0.5,
                                     z=vector[2]+0.5,
                                     text=f"{self.abc[index]}",
                                     showarrow=False))

        return go.Scatter3d(x=[0, vector[0]],
                            y=[0, vector[1]],
                            z=[0, vector[2]],
                            mode='lines',
                            name='cell',
                            line=dict(color='grey', width=5), showlegend=False
                            )

    def _cell_vectors_2d(self, vector, index):
        """Internal method to handle 3d unitcell vectors.

        Args:
            vector (list): end point for radius vector
            index (int): 0-a,1-b,2-c

        Returns:
            plotly.go: line and name of the vector
        """
        projected = self.apply_projection(vector)

        self.fig.add_annotation(x=projected[0]+0.2,
                                y=projected[1]+0.2,
                                text=f"{self.abc[index]}",
                                showarrow=False)

        return go.Scatter(
            x=[0, projected[0]],
            y=[0, projected[1]],
            mode='lines',
            name='cell',
            line=dict(color='grey', width=1), showlegend=False
        )

    def _bond_line_3d(self, bond_pair):
        """Internal method to draw bonds in 3d.

        Args:
            bond_pair (dict): instance of the bond_pair dict with info about this bond

        Returns:
            plotly.go : line instance
        """
        return go.Scatter3d(
            x=[bond_pair['begin'][0], bond_pair['end'][0]],
            y=[bond_pair['begin'][1], bond_pair['end'][1]],
            z=[bond_pair['begin'][2], bond_pair['end'][2]],
            mode='lines',
            opacity=0.4,
            name=f"{bond_pair['atom1']}-{bond_pair['atom2']}",
            text=f"{np.round(bond_pair['bond_length'],2)}",
            hoverinfo="text+name",
            line=dict(color=bond_pair['color'], width=7), showlegend=False
        )

    def _bond_line_2d(self, bond_pair):
        """Internal method to draw bonds in to 2d.

        Args:
            bond_pair (dict): instance of the bond_pair dict with info about this bond

        Returns:
            plotly.go : line instance
        """
        begin = self.apply_projection(bond_pair['begin'])
        end = self.apply_projection(bond_pair['end'])
        return go.Scatter(x=[begin[0], (begin[0]+end[0])/2, end[0]], y=[begin[1], (begin[1]+end[1])/2, end[1]],
                          mode='lines', opacity=0.4,
                          name=f"{bond_pair['atom1']}-{bond_pair['atom2']}",
                          text=f"{np.round(bond_pair['bond_length'],2)}",
                          hoverinfo="name+text",
                          line=dict(color=bond_pair['color'], width=7),
                          showlegend=False
                          )

    def _ms(self, pos, radius):
        """local method to draw atomic spheres

        Args:
            pos (list): postion of the center
            radius (float): radius of the sphere

        Returns:
            (x,y,z): mesh for the sphere
        """
        u, v = np.mgrid[0:2*np.pi:self.param['resolution']*2j,  # type: ignore
                        0:np.pi:self.param['resolution']*1j]  # type: ignore
        x = radius * np.cos(u)*np.sin(v) + pos[0]
        y = radius * np.sin(u)*np.sin(v) + pos[1]
        z = radius * np.cos(v) + pos[2]
        return (x, y, z)

    def info(self):
        """
        Display information about the atomic structure including attributes and methods.

        Returns:
        - str: A formatted string containing details about the structure.
        """
        info_str = "Object information:\n"
        
        # List of attributes
        attributes = [attr for attr in dir(self) if not attr.startswith('_') and not callable(getattr(self, attr))]
        
        # List of methods
        methods = [method for method in dir(self) if callable(getattr(self, method)) and not method.startswith('_')]
        
        info_str += "Attributes:\n"
        for attr in attributes:
            info_str += f"  - {attr}: {getattr(self, attr)}\n"
        
        info_str += "Methods:\n"
        for method in methods:
            info_str += f"  - {method}\n"
        
        print(info_str)