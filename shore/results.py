import os
import re
import numpy as np
import pickle
import plotly.graph_objects as go

class ResultsHandler:
    def __init__(self, path=None, name=None, load=None):
        """Initialize with the path and filename."""
        self.name = name
        self.path = path
        if name and path:
            self.save_file_path=f'{self.path}/{self.name}.pkl'
        self.data = {}
        self.available = {
            "site_number": set(),
            "core_level": set(),
            "polarization": set(),
            "element": set()
        }
        if self.path:
            self.files = self._find_files()
            for file in self.files:
                self._parse_data(file)
        if load:
            self._load(file=load)

    def _find_files(self):
        """Find all files starting with 'absspct' in the specified folder."""
        files = []
        try:
            for filename in os.listdir(self.path):
                if filename.startswith("absspct"):
                    files.append(filename)
        except FileNotFoundError:
            print(f"Error: The directory '{self.path}' does not exist.")
        return files

    def _parse_data(self, file):
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

    def get_data(self, element=None, edge=None, site=None, polarization=None):
        """Retrieve data based on specified parameters."""
        try:
            return self.data[element][edge][site][polarization]
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
            
  

    def plot(self, fig=None, element=None, core_level=None, site_number=None, polarization=None, name=None, norm=True, lw=2, lc=None):
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
