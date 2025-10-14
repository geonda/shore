import os
import re
import numpy as np
import pickle
import plotly.graph_objects as go
from plotly.subplots import make_subplots

class ResultsHandler:
    def __init__(self, path=None, name=None, load=None):
        """Initialize with the path and filename."""
        self.name = name
        self.path = path
        if name and path:
            self.save_file_path=f'{self.path}/{self.name}.pkl'
        self.data = {}
        self.data_rixs={}
        self.available = {
            "site_number": set(),
            "core_level": set(),
            "polarization": set(),
            "element": set(),
            'photon_in':set(),
            'photon_out':set(),
            'energy_point':set(),
        }
        if self.path:
            self.files_xas,self.files_rxs = self._find_files()
            for file in self.files_xas:
                self._parse_data_xas(file)
            for file in self.files_rxs:
                self._parse_data_rxs(file)
        if load:
            self._load(file=load)

    def _find_files(self):
        """Find all files starting with 'absspct' in the specified folder."""
        files_xas = []
        files_rxs = []
        try:
            for filename in os.listdir(self.path):
                if filename.startswith("absspct"):
                    files_xas.append(filename)
            for filename in os.listdir(self.path):
                if filename.startswith("rxsspct"):
                    files_rxs.append(filename)
        except FileNotFoundError:
            print(f"Error: The directory '{self.path}' does not exist.")
        return files_xas, files_rxs
    
    def _parse_data_rxs(self,file):
        
        pattern = r"rxsspct_(\w+)[._]?(\w+)_(\d+)\.(\d+)\.(\d+)"

        match = re.match(pattern, file)

        if match:
            element = match.group(1)        # F
            core_level = match.group(2)     # 1s
            photon_in = int(match.group(3))      # 01
            energy_point = int(match.group(4))  # 0001 -> 1
            photon_out = int(match.group(5))     # 02

            # Update available parameters
            self.available["element"].add(element)
            self.available["core_level"].add(core_level)
            self.available["photon_in"].add(photon_in)
            self.available["energy_point"].add(energy_point)
            self.available["photon_out"].add(photon_out)

            # Initialize nested dictionaries if needed
            if element not in self.data_rixs:
                self.data_rixs[element] = {}
            if core_level not in self.data_rixs[element]:
                self.data_rixs[element][core_level] = {}
            if photon_in not in self.data_rixs[element][core_level]:
                self.data_rixs[element][core_level][photon_in] = {}
            if energy_point not in self.data_rixs[element][core_level][photon_in]:
                self.data_rixs[element][core_level][photon_in][energy_point] = {}
            if photon_out not in self.data_rixs[element][core_level][photon_in][energy_point]:
                self.data_rixs[element][core_level][photon_in][energy_point][photon_out] = {}

            # Load data from file and store
            tmp = np.loadtxt(os.path.join(self.path, file)).transpose()
            self.data_rixs[element][core_level][photon_in][energy_point][photon_out]['energy'] = tmp[0]
            self.data_rixs[element][core_level][photon_in][energy_point][photon_out]['spectrum'] = tmp[2]

        else:
            raise ValueError(f"Filename '{file}' does not match expected rxsspct")

    def _parse_data_xas(self, file):
        """Parse the filename to extract element, site number, core-level, and polarization."""
        # Regex pattern to match the expected filename format
        pattern = r"absspct_(\w+)\.(\d+)_(\w+)_(\d+)"
        match = re.match(pattern, file)

        if match:
            element = match.group(1)  # Extract element (e.g., Ti)
            site_number = int(match.group(2))  # Extract site number (e.g., 0001)
            core_level = match.group(3)  # Extract core level (e.g., 2p)
            polarization = int(match.group(4))  # Extract polarization (e.g., 02)

            # Update available parameters
            self.available["element"].add(element)
            self.available["site_number"].add(site_number)
            self.available["core_level"].add(core_level)
            self.available["polarization"].add(polarization)

            # Initialize nested dictionaries if they don't exist
            if element not in self.data:
                self.data[element] = {}
            if core_level not in self.data[element]:
                self.data[element][core_level] = {}
            if site_number not in self.data[element][core_level]:
                self.data[element][core_level][site_number] = {}

            # Load data from the file and store it in the nested dictionary
            tmp = np.loadtxt(os.path.join(self.path, file)).transpose()
            if polarization not in self.data[element][core_level][site_number]:
                self.data[element][core_level][site_number][polarization] = {}
                
            self.data[element][core_level][site_number][polarization]['energy'] = tmp[0]
            self.data[element][core_level][site_number][polarization]['spectrum'] = tmp[1]
            
        else:
            raise ValueError(f"Filename '{file}' does not match expected format.")

    def get_data_xas(self, element=None, edge=None, site=None, polarization=None):
        """Retrieve data based on specified parameters."""
        try:
            return self.data[element][edge][site][polarization]
        except KeyError as e:
            print(f"Error: {e} - Please check your parameters.")
    
    def get_data_rixs(self, element=None, edge=None,  photon_in=None, energy_point=None, photon_out=None):
        """Retrieve data based on specified parameters."""
        try:
            return self.data_rixs[element][edge][photon_in][energy_point][photon_out]
        except KeyError as e:
            print(f"Error: {e} - Please check your parameters.")

    def save(self):
        """Save the current state of the object to a pickle file."""
        with open(os.path.join(self.path, f'{self.name}.pkl'), 'wb') as f:
            pickle.dump(self, f)
    
    def _load(self, file=None):
        """
        Load the object state from a pickle file.

        Args:
            file (str): Optional; path to the pickle file. If not provided, 
                         it will load from '{self.path}/{self.name}.pkl'.

        Returns:
            None: The method populates the object's attributes with 
                  the loaded data.
        """
        # Determine the file path
        if not file:
            file = os.path.join(self.path, f'{self.name}.pkl')

        # Load the object state from the specified pickle file
        with open(file, 'rb') as f:
            loaded_instance = pickle.load(f)

            # Populate current instance's attributes with loaded instance's attributes
            for attr in vars(loaded_instance):
                setattr(self, attr, getattr(loaded_instance, attr))

    def load(self, file=None):
        if not file:
            """Load the object state from a pickle file."""
            with open(os.path.join(self.path, f'{self.name}.pkl'), 'rb') as f:
                loaded_instance = pickle.load(f)
                return loaded_instance
        else:
            """Load the object state from a pickle file."""
            with open(f'{file}', 'rb') as f:
                loaded_instance = pickle.load(f)
                return loaded_instance

    def plot_rixs(self, fig=None, element=None, core_level=None, photon_in=None, photon_out=None, energy_point=None, name=None, norm=True, lw=2, lc=None):
        """Plot RIXS vs energy loss using Plotly."""
        if not fig:
            fig = go.Figure()

        # Check if specific parameters are provided
        if element is None:
            elements = self.data.keys()  # Get all available elements
        else:
            elements = [element]

        if core_level is None:
            core_levels = set()
            for el in elements:
                core_levels.update(self.data[el].keys())  # Collect all available core levels
        else:
            core_levels = [core_level]
        if not photon_in:
            photon_in=list(self.available['photon_in'])
        else:
            photon_in=list(photon_in)

        if not photon_out:
            photon_out=list(self.available['photon_out'])
        else:
            photon_out=list(photon_out)
        if not energy_point:
            energy_point=1

        # Loop through the selected elements, core levels, and site numbers
        for el in elements:
            for cl in core_levels:
                spectrum=0
                for photon_ini in photon_in:
                    for photon_outi in photon_out:
                        spectrum += self.data_rixs[el][cl][photon_ini][energy_point][photon_outi]['spectrum']
                        energy=self.data_rixs[el][cl][photon_ini][energy_point][photon_outi]['energy']
                if norm:  # Normalize the spectrum if required
                    max_intensity = max(spectrum)  # Avoid division by zero
                    spectrum = [s / max_intensity for s in spectrum]

                if not name:
                    name = f'{el} Energy in point: {energy_point}'  # Updated name without polarization
                    
                if lc:
                    # Add trace for each polarization
                    fig.add_trace(go.Scatter(
                        x=energy,
                        y=spectrum,
                        mode='lines',
                        name=name,
                        line=dict(width=lw,color=lc)
                    ))
                else:
                    fig.add_trace(go.Scatter(
                        x=energy,
                        y=spectrum,
                        mode='lines',
                        name=name,
                        line=dict(width=lw)
            ))
                    
    

        # Update layout for scientific styling with legend settings
        fig.update_layout(
            # title='X-ray Absorption Spectrum',
            xaxis_title='Energy  Loss (eV)',
            yaxis_title='RIXS Intensity (arb. units)',
            template='plotly_white',
            font=dict(size=12),
            hovermode='x unified',
            showlegend=True  # Ensure legend is shown
        )
        from shore import plotly_formatting
        fig=plotly_formatting(fig)
        fig.update_xaxes(range=[0,30])
        return fig

    def readout_rixsmap(self,element='O', core_level='1s', photon_in=1, photon_out=2, starting_energy=-2, delta_energy=2):
        z = []
        y = []
        for step in sorted(list(self.data_rixs[element][core_level][photon_in].keys())):
            z.append(self.data_rixs[element][core_level][photon_in][step][photon_out]['spectrum'])
            y.append(starting_energy+delta_energy*(float(step)-1))
        x = self.data_rixs[element][core_level][photon_in][step][photon_out]['energy']
        return x, y, z
            
    def rixsmap(self,xas=None,element='O', core_level='1s', photon_in=1, photon_out=2, starting_energy=-2, delta_energy=2, rshift=0, ):


        custom_colorscale=[
            [0, 'blue'],  # Start color (blue)
            [0.5, 'white'],  # Middle color (white)
            [1, 'red']  # End color (red)
            ]
        x,y,z=self.readout_rixsmap(element=element, core_level=core_level, 
                                    photon_in=photon_in, 
                                    photon_out=photon_out, 
                                    starting_energy=starting_energy,
                                     delta_energy=delta_energy)
        y=[yi+rshift for yi in y]
        contour_fig = go.Figure(go.Heatmap(
            x=x, y=y, z=z, colorscale=custom_colorscale, colorbar={"title": "RIXS"}, zsmooth='best', opacity=0.6))
        if xas:
            xas_spectrum=0
            for site in sorted(list(xas.data[element][core_level].keys())):
                xas_spectrum+=xas.data[element][core_level][site][int(photon_in)]['spectrum']
            xas_energy = xas.data[element][core_level][site][int(photon_in)]['energy']
            # Create the scatter plot
            scatter_fig = go.Figure(go.Scatter(
                x=xas_spectrum, y=xas_energy+rshift, mode='lines', ))
        else:
            xas_spectrum=0
            for site in sorted(list(self.data[element][core_level].keys())):
                xas_spectrum+=self.data[element][core_level][site][int(photon_in)]['spectrum']
            xas_energy = self.data[element][core_level][site][int(photon_in)]['energy']
            # Create the scatter plot
            scatter_fig = go.Figure(go.Scatter(
                x=xas_spectrum, y=xas_energy+rshift, mode='lines', ))

        # Create subplots with shared x-axis
        fig = make_subplots(rows=1, cols=2, shared_yaxes=True,
                            column_widths=[0.2, 0.8], horizontal_spacing=0.01)

        # Add contour plot to the first subplot
        for trace in contour_fig.data:
            fig.add_trace(trace, row=1, col=2)

        # Add scatter plot to the second subplot
        for trace in scatter_fig.data:
            fig.add_trace(trace, row=1, col=1)

        fig.update_layout(xaxis_title='XAS', yaxis_title='Energy, eV',
                        margin=dict(l=0, r=0, b=0, t=1))
        fig.update_layout(xaxis2_title='Energy loss, eV', yaxis_title='Energy in, eV',
                        margin=dict(l=0, r=0, b=0, t=1))
        fig.update_xaxes(autorange="reversed", row=1, col=1)
        fig.update_layout(xaxis2=dict(range=[0, 20]))
        fig.update_layout(yaxis2=dict(range=[min(y), max(y)]))
        fig.update_layout(yaxis1=dict(range=[min(y), max(y)]))
        # Set mirrored axes
        fig.update_layout(xaxis=dict(mirror=True, color='black', showline=True,
                                    linewidth=2, linecolor='black', showticklabels=False), yaxis=dict(mirror=True, showline=True, linewidth=2, linecolor='black'))
        fig.update_layout(xaxis2=dict(mirror=True, color='black', showline=True,
                                    linewidth=2, linecolor='black'), yaxis2=dict(mirror=True, showline=True, linewidth=2, linecolor='black'))
        fig.update_xaxes(zeroline=True, zerolinewidth=1, zerolinecolor='black')

        # fig.add_hline(y=529.7, line_dash="dot", row=1, col=1,
        #             line_color="black", line_width=2)
        # fig.add_hline(y=529.7, line_dash="dot", row=1, col=2,
        #             line_color="black", line_width=2)
        # Set background color to white
        fig.update_layout(plot_bgcolor='white', paper_bgcolor='white', font=dict(
            size=20,
        ))
        return fig


    def plot_xas(self, fig=None, element=None, core_level=None, site_number=None, polarization=None, name=None, norm=True, lw=2, lc=None):
        """Plot XAS absorption vs energy using Plotly."""
        if not fig:
            fig = go.Figure()

        # Check if specific parameters are provided
        if element is None:
            elements = self.data.keys()  # Get all available elements
        else:
            elements = [element]

        if core_level is None:
            core_levels = set()
            for el in elements:
                core_levels.update(self.data[el].keys())  # Collect all available core levels
        else:
            core_levels = [core_level]

        if site_number is None:
            site_numbers = set()
            for el in elements:
                for cl in core_levels:
                    site_numbers.update(self.data[el][cl].keys())  # Collect all available site numbers
        else:
            site_numbers = [site_number]

        # Loop through the selected elements, core levels, and site numbers
        for el in elements:
            for cl in core_levels:
                for sn in site_numbers:
                    spectrum = 0
                    if polarization is None:
                        # Sum over all polarizations
                        for pol in self.data[el][cl][sn]:
                            energy = self.data[el][cl][sn][pol]['energy']
                            spectrum += self.data[el][cl][sn][pol]['spectrum']
                        
                        if norm:  # Normalize the spectrum if required
                            max_intensity = max(spectrum)  # Avoid division by zero
                            spectrum = [s / max_intensity for s in spectrum]
                        
                        if not name:
                            name = f'{el} Site {sn}'  # Updated name without polarization
                        
                        if lc:
                            # Add trace for each polarization
                            fig.add_trace(go.Scatter(
                                x=energy,
                                y=spectrum,
                                mode='lines',
                                name=name,
                                line=dict(width=lw,color=lc)
                            ))
                        else:
                            fig.add_trace(go.Scatter(
                                x=energy,
                                y=spectrum,
                                mode='lines',
                                name=name,
                                line=dict(width=lw)
                            ))
                    else:
                        # Specific polarization case
                        energy = self.data[el][cl][sn][polarization]['energy']
                        spectrum += self.data[el][cl][sn][polarization]['spectrum']

                        if norm:  # Normalize the spectrum if required
                            max_intensity = max(spectrum)  # Avoid division by zero
                            spectrum = [s / max_intensity for s in spectrum]
                        
                        if not name:
                            name = f'{el} {cl} Site {sn} Polarization {polarization}'
                        
                        # Add trace for the specified polarization
                        if lc:
                            # Add trace for each polarization
                            fig.add_trace(go.Scatter(
                                x=energy,
                                y=spectrum,
                                mode='lines',
                                name=name,
                                line=dict(width=lw,color=lc)
                            ))
                        else:
                            fig.add_trace(go.Scatter(
                                x=energy,
                                y=spectrum,
                                mode='lines',
                                name=name,
                                line=dict(width=lw)
                            ))

        # Update layout for scientific styling with legend settings
        fig.update_layout(
            # title='X-ray Absorption Spectrum',
            xaxis_title='Relative Energy (eV)',
            yaxis_title='XAS Intensity (arb. units)',
            template='plotly_white',
            font=dict(size=12),
            hovermode='x unified',
            showlegend=True  # Ensure legend is shown
        )
        from shore import plotly_formatting
        fig=plotly_formatting(fig)
        return fig

                # Show the figure
                

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
