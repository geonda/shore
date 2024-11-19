
import json

import numpy as np

from ase.atoms import Atoms
  # Assuming you may need this function later
from ase.units import Bohr


import spglib

from IPython.display import JSON as js

# Shore-related imports


# Optional: If you plan to use threading, you can keep this import.
import threading


{

        "dft.program": "qe",
        'dft.den.kmesh':  '-2',
        "dft.ecut":' -1',

        "screen.nbands":  10,
        "screen.final.dr": "0.02",
        "screen.shells": "6.0",
        "screen.kmesh": '-2',
        "screen.core_offset.enable" : "true",

        'bse.core.broaden': '-1',
        'bse.core.haydock.converge.thresh' : '0.001 5',
        'bse.core.screen_radius': ' 5.5',
        'bse.nbands':  10,
        "bse.core.haydock.niter":  '1000',
        "bse.kmesh":' -2',

        "psp.pp_database": "ONCVPSP-PBE-PDv0.4-stringent",
        "psp.ecut_quality":  'high',
        "opf.program": "hamann",
        "computer.para_prefix": 'srun ',

        }


class default:
    def __init__(self):
        self.input={


    #dft
    "dft.program": "qe",
    'dft.den.kmesh':  '-2',
    "dft.ecut":' -1',
    "diemac":10000,

    # screen
    "screen.nbands":  10,
    "screen.kmesh": '-2',

    "screen.shells": "6.0",
    "screen.core_offset.enable" : "true",
    "screen.final.dr": "0.02",
    "screen.grid.rmax": "10",
    "screen.grid.rmode": "lagrange uniform",
    "screen.grid.ang": "5 11 11 9 7",
    "screen.grid.deltar": "0.10 0.15 0.25 0.25 0.25",
    "screen.grid.shells": " -1 4 6 8 10",
    "screen.lmax": "2",
    "screen.shells": "3.5 4.0 4.5 5.0 5.5 6.0",

    #bse
    "bse.kmesh":' -2',
    'bse.nbands':  10,

    'bse.core.broaden': '-1',
    'bse.core.haydock.converge.thresh' : '0.001 5',
    'bse.core.screen_radius': ' 5.5',
    "bse.core.haydock.niter":  '1000',
    "cnbse.spect_range": "1000 -10 40",

    "psp.pp_database": "ONCVPSP-PBE-PDv0.4-stringent",
    "psp.ecut_quality":  'high',
    "opf.program": "hamann",
    "computer.para_prefix": 'mpirun -np 4 ',

}




class EnergyCalculator:
    def __init__(self, element=None, edge=None):
        self.edge=edge
        self.element=element

        # Dictionary to hold edge-energy pairs
        self.edge_energy={}
        self.edge_energy["K"] = {
                            "H": 13.6,
                            "He": 24.6,
                            "Li": 54.7,
                            "Be": 111.5,
                            "B": 188,
                            "C": 284.2,
                            "N": 409.9,
                            "O": 543.1,
                            "F": 696.7,
                            "Ne": 870.2,
                            "Na": 1070.8,
                            "Mg": 1303,
                            "Al": 1559,
                            "Si": 1839,
                            "P": 2145.5,
                            "S": 2472,
                            "Cl": 2822.4,
                            "Ar": 3205.9,
                            "K": 3608.4,
                            "Ca": 4038.5,
                            "Sc": 4492,
                            "Ti": 4966,
                            "V": 5465,
                            "Cr": 5989,
                            "Mn": 6539,
                            "Fe": 7112,
                            "Co": 7709,
                            "Ni": 8333,
                            "Cu": 8979,
                            "Zn": 9659,
                            "Ga": 10367,
                            "Ge": 11103,
                            "As": 11867,
                            "Se": 12658,
                            "Br": 13474,
                            "Kr": 14326,
                            "Rb": 15200,
                            "Sr": 16105,
                            "Y": 17038,
                            "Zr": 17998,
                            "Nb": 18986,
                            "Mo": 20000,
                            "Tc": 21044,
                            "Ru": 22117,
                            "Rh": 23220,
                            "Pd": 24350,
                            "Ag": 25514,
                            "Cd":26711, 
                            "In" :27940, 
                            "Sn" :29200, 
                            "Sb" :30491, 
                            "Te" :31814, 
                            "I" :33169, 
                            "Xe" :34561, 
                            "Cs" :35985, 
                            "Ba" :37441, 
                            "La" :38925, 
                            "Ce" :40443, 
                            "Pr" :41991, 
                            "Nd" :43569, 
                            "Pm" :45184, 
                            "Sm" :46834, 
                            "Eu" :48519, 
                            "Gd" :50239, 
                            "Tb" :51996, 
                            "Dy" :53789, 
                            "Ho" :55618, 
                            "Er" :57486, 
                            "Tm" :59390, 
                            "Yb" :61332, 
                            "Lu" :63314, 
                            "Hf" :65351, 
                            "Ta" :67416, 
                            "W" :69525, 
                            "Re" :71676, 
                            "Os" :73871, 
                            "Ir" :76111, 
                            "Pt" :78395, 
                            "Au" :80725, 
                            "Hg" :83102, 
                            "Tl" :85530, 
                            "Pb ":88005,  
                            "Bi ":90526,  
                            "Po ":93105,  
                            "At ":95730,  
                            "Rn ":98404,  
                            "Fr":101137,  
                            "Ra":103922,  
                            "Ac":106755,  
                            "Th":109651,  
                            "Pa":112601,  
                            "U":115606}
        self.edge_energy["L1"]={
                            "N": 37.3,
                            "O": 41.6,
                            "F": 0.0,
                            "Ne": 48.5,
                            "Na": 63.5,
                            "Mg": 88.6,
                            "Al": 117.8,
                            "Si": 149.7,
                            "P": 189,
                            "S": 230.9,
                            "Cl": 270.0,
                            "Ar": 326.3,
                            "K": 378.6,
                            "Ca": 438.4,
                            "Sc": 498,
                            "Ti": 560.9,
                            "V": 626.7,
                            "Cr": 696,
                            "Mn": 769.1,
                            "Fe": 844.6,
                            "Co": 925.1,
                            "Ni": 1008.6,
                            "Cu": 1096.7,
                            "Zn": 1196.2,
                            "Ga": 1299,
                            "Ge": 1414.6,
                            "As": 1527,
                            "Se": 1652,
                            "Br": 1782,
                            "Kr": 1921,
                            "Rb": 2065,
                            "Sr": 2216,
                            "Y": 2373,
                            "Zr": 2532,
                            "Nb": 2698,
                            "Mo": 2866,
                            "Tc": 3043,
                            "Ru": 3224,
                            "Rh": 3412,
                            "Pd": 3604,
                            "Ag": 3806,
                            "Cd": 4018,
                            "In":4238, 
                            "Sn" :4465, 
                            "Sb" :4698, 
                            "Te" :4939, 
                            "I" :5188, 
                            "Xe" :5453, 
                            "Cs" :5714, 
                            "Ba" :5989, 
                            "La" :6266, 
                            "Ce" :6548, 
                            "Pr" :6835, 
                            "Nd" :7126, 
                            "Pm" :7428, 
                            "Sm" :7737, 
                            "Eu" :8052, 
                            "Gd" :8376, 
                            "Tb ":8708, 
                            "Dy ":9046, 
                            "Ho ":9394, 
                            "Er ":9751, 
                            "Tm ":10116, 
                            "Yb ":10486, 
                            "Lu ":10870, 
                            "Hf ":11271, 
                            "Ta ":11682, 
                            "W ":12100, 
                            "Re ":12527, 
                            "Os ":12968, 
                            "Ir ":13419, 
                            "Pt ":13880, 
                            "Au ":14353, 
                            "Hg ":14839, 
                            "Tl ":15347, 
                            "Pb ":15861, 
                            "Bi ":16388, 
                            "Po ":16939, 
                            "At ":17493, 
                            "Rn ":18049, 
                            "F r ":18639, 
                            "Ra ":19237 ,  
                            "Ac":19840 ,  
                            "Th":20472 ,  
                            "Pa":21105 ,  
                            "U" :21757
                            }
        self.edge_energy["L2"]={
                        "F": 48.5,
                        "Ne": 63.5,
                        "Na": 88.6,
                        "Mg": 17.8,
                        "Al": 49.7,
                        "Si": 89,
                        "P": 136,
                        "S": 163.6,
                        "Cl": 202,
                        "Ar": 250.6,
                        "K": 297.3,
                        "Ca": 349.7,
                        "Sc": 403.6,
                        "Ti": 460.2,
                        "V": 519.8,
                        "Cr": 583.8,
                        "Mn": 649.9,
                        "Fe": 719.9,
                        "Co": 793.2,
                        "Ni": 870,
                        "Cu": 952.3,
                        "Zn": 1044.9,
                        "Ga": 1143.2,
                        "Ge": 1248.1,
                        "As": 1359.1,
                        "Se": 1474.3,
                        "Br": 1596,
                        "Kr": 1730.9,
                        "Rb": 1864,
                        "Sr": 2007,
                        "Y": 2156,
                        "Zr": 2307,
                        "Nb": 2465,
                        "Mo": 2625,
                        "Tc":2793, 
                        "Ru" :2967, 
                        "Rh" :3146, 
                        "Pd" :3330, 
                        "Ag" :3524, 
                        "Cd" :3727, 
                        "In" :3938, 
                        "Sn" :4156, 
                        "Sb" :4380, 
                        "Te" :4612, 
                        "I" :4852, 
                        "Xe" :5107, 
                        "Cs" :5359, 
                        "Ba" :5624, 
                        "La" :5891, 
                        "Ce" :6164, 
                        "Pr" :6440, 
                        "Nd" :6722, 
                        "Pm" :7013, 
                        "Sm" :7312, 
                        "Eu" :7617, 
                        "Gd" :7930, 
                        "Tb ":8252, 
                        "Dy ":8581, 
                        "Ho ":8918, 
                        "Er ":9264, 
                        "Tm ":9617, 
                        "Yb ":9978, 
                        "Lu ":10349, 
                        "Hf ":10739, 
                        "Ta ":11136, 
                        "W ":11544 ,  
                        "Re":11959 ,  
                        "Os":12385 ,  
                        "Ir":12824 ,  
                        "Pt":13273 ,  
                        "Au":13734 ,  
                        "Hg":14209 ,  
                        "Tl":14698 ,  
                        "Pb":15200 ,  
                        "Bi":15711 ,  
                        "Po":16244 ,  
                        "At":16785 ,  
                        "Rn":17337 ,  
                        "Fr":17907 ,  
                        "Ra":18484 ,  
                        "Ac":19083 ,  
                        "Th":19693 ,  
                        "Pa":20314 ,  
                        "U" :20948
                        }



        
    def get_energy(self):
        """Returns the energy corresponding to the specified edge."""
        return self.edge_energy[self.edge].get(self.element, "Edge not found")

class Photon:
    def __init__(self,interaction='dipole', polarization=[0, 0, 0], q=[0, 0, 1], energy=None):
        self.polarization = polarization
        self.q = q
        self.element=energy['element']
        self.edge=energy['edge']
        self.photon_energy = EnergyCalculator(edge=self.edge,element=self.element).get_energy()
        self.interaction=interaction

    def info(self):
        """Return a string with information about the photon."""
        return f"Photon(polarization={self.polarization}, q={self.q}, energy={self.photon_energy} eV, light-matter interaction={self.interaction})"

class Light:
    def __init__(self, photons=[]):
        if not photons:
            self.photons = []
        else:
            self.photons=photons
        

    def add(self, photon):
        """Add a photon to the light bundle."""
        if isinstance(photon, Photon):
            self.photons.append(photon)
        else:
            raise ValueError("Only instances of Photon can be added.")

    def show(self):
        """Display information about all photons in the light bundle and total energy."""
        if not self.photons:
            print("The light bundle is empty.")
            return
        
        
        print("Light Bundle Information:")
        for idx, photon in enumerate(self.photons):
            print(f"{idx + 1}: {photon.info()}")

    def write_to_folder(self, file_path='./', ):
        """Write each photon's information to a separate file."""
        for idx, photon in enumerate(self.photons):
            filename = f"photon{idx + 1}"
            with open(f"{file_path}/{filename}", 'w') as f:
                f.write(f"{photon.interaction}\n")
                f.write(f"cartesian  {photon.polarization[0]} {photon.polarization[1]} {photon.polarization[2]}\n")
                f.write("end\n")
                f.write(f"cartesian  {photon.q[0]} {photon.q[1]} {photon.q[2]}\n")
                f.write("end\n")
                f.write(f"{photon.photon_energy}\n")

class Matter:
    def __init__(self,structure=None, **kwargs):
        """
        Initializes the Matter class with DFT, CNBSE, and general dictionaries.
       
        """
        self.structure=structure
        if 'load' in kwargs.keys():
            self.params=self.read_from_file(kwargs['load'])
            del kwargs['load']
        else:
            self.params=default().input
        for k,v in kwargs.items():
            if k=='screen_dft_energy_range':
                k='screen.dft_energy_range'
            elif  k=='bse_dft_energy_range':
                k='bse.dft_energy_range'
            elif k=="computer_para_prefix":
                k="computer.para_prefix"
            if '_' in k:
                k=k.replace('_','.')
            self.params.update({k:v})
        # print(self.params)
        for item in self.params.keys():
            setattr(self,item,self.params[item])
        

        
    def show(self):
        return js(self.params)

    def read_from_file(self, file_path):
        tmp = {}
        
        try:
            with open(file_path, 'r') as file:
                lines = file.readlines()
                global_index = 0
                
                while global_index < len(lines):
                    # Strip whitespace and ignore empty lines
                    line = lines[global_index].strip()
                    # print(line)
                    # if '#' in line:
                    #     global_index += 1
                    #     continue
                    if not line:
                        global_index += 1  # Increment index to avoid infinite loop
                        continue
                    
                    # print(line)  # For debugging
                    
                    # Check if both '{' and '}' are in the line
                    if '{' in line and '}' in line:
                        key, value = line.split('{', 1)
                        key = key.strip()
                        value = value.strip(' }')  # Remove closing brace and whitespace
                        tmp[key] = value
                    
                    else: 
                        key = line.split('{')[0].strip()
                        value = []
                        local_index = 0
                        
                        while True:
                            cline = lines[global_index + local_index].strip()
                            if '{' in cline:
                                if cline.split('{')[1].strip():
                                    value.append(cline.split('{')[1].strip())
                            elif '}' not in cline:
                                value.append(cline)
                            else:
                                # Split on the first occurrence of '}' to capture the last part
                                if cline.split('}')[0].strip():
                                    value.append(cline.split('}')[0].strip())
                                break  # Exit loop when closing brace is found
                            local_index += 1
                        
                        global_index += local_index  # Move the global index forward
                        tmp[key] = value  # Join multiline values into a single string

                    global_index += 1  # Move to the next line
        except Exception as e:
            print(f"An error occurred: {e}")
        return tmp



        
    def info(self):
        """
        Display information about the atomic structure including attributes and methods.

        Returns:
        - str: A formatted string containing details about the structure.
        """
        info_str = "Atomic Structure Information:\n"
        
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


class Input():
    def __init__(self, target='xas', name='ocean', matter=None, light=None):
        self.name=name
        self.target=target
        self.structure=matter.structure
        self.matter=matter
        self.light=light

        if 'xas' in target:
            self.content=ocean_input(matter=self.matter, light=self.light, )
        elif 'qe' in target:
            print('qe_input')
            # self.data=qe_input()
        else:
            print('Wrong target')

    def fork(self,**kwargs):
        import copy
        tmp=copy.deepcopy(self)
        for k,v in kwargs.items():
            if k=='name':
                tmp.name=v
            elif k=='screen_dft_energy_range':
                k='screen.dft_energy_range'
                tmp.content.input.update({k:v})
            elif  k=='bse_dft_energy_range':
                k='bse.dft_energy_range'
                tmp.content.input.update({k:v})
            elif '_' in k:
                k=k.replace('_','.')
                tmp.content.input.update({k:v})
            else:
                tmp.content.input.update({k:v})
        return tmp

    def info(self):
        """
        Display information about the atomic structure including attributes and methods.

        Returns:
        - str: A formatted string containing details about the structure.
        """
        info_str = "Atomic Structure Information:\n"
        
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
    
class ocean_input():
    def __init__(self, matter=None,  light=None):
        self.structure=matter.structure
        if self.structure:
            self.atoms=self.structure.ocean_atoms
            self.structure_formatted=self._structure2input()
        self.matter=matter
        self.light=light
        self.edge=None
        self.element=None
        if light:
            self.edge_line=self._get_edge_line()
        try:
            self.input={**self.matter.params,**self.edge_line, **self.structure_formatted}
        except Exception as e:
            print("Empty input")

    def _get_types(self,species):
        types = []
        for Z in self.atoms.numbers:
            for n, Zs in enumerate(species):
                if Z == Zs:
                    types.append(n + 1)
        return types
    
    def _get_atomic_positions(self):
        atomic_positions_str = []
        for atom in self.atoms:
            atomic_positions_str.append( '{coords[0]:.10f} {coords[1]:.10f} {coords[2]:.10f}\n'.format(
                coords=[atom.a, atom.b, atom.c] ))
        return atomic_positions_str
        
    def _structure2input(self):
        species = sorted(set(self.atoms.numbers))
        structure_formatted=dict(
        znucl=' '.join(str(Z) for Z in species),
        typat=' '.join(str(item) for item in self._get_types(species)),
        xred=''.join(str(item) for item in self._get_atomic_positions()),
        acell='{acell[0]} {acell[0]} {acell[0]}'.format( acell=[1/Bohr] ) ,
        rprim=' {cell[0][0]:.14f} {cell[0][1]:.14f} {cell[0][2]:.14f}\n {cell[1][0]:.14f} {cell[1][1]:.14f} {cell[1][2]:.14f}\n {cell[2][0]:.14f} {cell[2][1]:.14f} {cell[2][2]:.14f} '.format(cell=self.atoms.cell),
        )
        return structure_formatted

    def _get_edge_line(self):
        try: 
            photon=self.light.photons[0]
            self.element=photon.element
            self.edge = photon.edge
            edge_line=f"-{self._get_element_number(self.element)} {self._get_edge_number(self.edge)}"
        except:
            print(
                'something wrong with light'
            )
            edge_line=None
        
        return dict(edges=edge_line) 

    def _get_element_number(self, element):
        for s,n  in zip(self.structure.atoms.get_chemical_symbols(),self.structure.atoms.get_atomic_numbers()):
            if s==element:
                return n
        print('Error wrong elements perhaps')
        return None

    def _get_edge_number(self,edge):
        if 'K' in edge:
            return '1 0'
        elif 'L' in edge:
            return '2 1'
        else: 
            print("wrong edge perhaps")
            return None

    def read_from_file(self, file_path):
        tmp = {}
        
        try:
            with open(file_path, 'r') as file:
                lines = file.readlines()
                global_index = 0
                
                while global_index < len(lines):
                    # Strip whitespace and ignore empty lines
                    line = lines[global_index].strip()
                    print(line)
                    # if '#' in line:
                    #     global_index += 1
                    #     continue
                    if not line:
                        global_index += 1  # Increment index to avoid infinite loop
                        continue
                    
                    # print(line)  # For debugging
                    
                    # Check if both '{' and '}' are in the line
                    if '{' in line and '}' in line:
                        key, value = line.split('{', 1)
                        key = key.strip()
                        value = value.strip(' }')  # Remove closing brace and whitespace
                        tmp[key] = value
                    
                    else: 
                        key = line.split('{')[0].strip()
                        value = []
                        local_index = 0
                        
                        while True:
                            cline = lines[global_index + local_index].strip()
                            if '{' in cline:
                                if cline.split('{')[1].strip():
                                    value.append(cline.split('{')[1].strip())
                            elif '}' not in cline:
                                value.append(cline)
                            else:
                                # Split on the first occurrence of '}' to capture the last part
                                if cline.split('}')[0].strip():
                                    value.append(cline.split('}')[0].strip())
                                break  # Exit loop when closing brace is found
                            local_index += 1
                        
                        global_index += local_index  # Move the global index forward
                        tmp[key] = value  # Join multiline values into a single string

                    global_index += 1  # Move to the next line
        except Exception as e:
            print(f"An error occurred: {e}")
        self.input=tmp
        try:
            rprim=np.array([[float(f) for f in tmp['rprim'][j].split()] for j in range(3)])
            xred=np.array([[float(f) for f in tmp['xred'][j].split()] for j in range(len(tmp['xred']))])
            typat=tmp['typat'].split()
            znucl=tmp['znucl'].split()
            positions = []
            atoms_=[int(znucl[int(i)-1]) for i in typat]
            print(atoms_)
            # Convert fractional coordinates to Cartesian coordinates
            for i, frac_coord in enumerate(xred):
                # Calculate Cartesian coordinates based on rprim and acell
                cartesian_coord = np.dot(frac_coord, rprim) * np.array([float(f) for f in tmp['acell'][0].split()])
                positions.append(cartesian_coord)
            print(positions)
            atoms = Atoms(symbols=atoms_, positions=positions)
            self.structure=atoms
            # self.matter.params={key: tmp[key] for key in tmp.keys() if key not in ['rprim', 'typat', 'acell', 'xred', 'znucl']}
        except Exception as e:
            print(e)
        
    def show(self):
        return js(self.input)

    def _flatten_dict(self,d):
        items = {}
        for k, v in d.items():
            new_key = f"{k}"
            if isinstance(v, dict):
                items.update(self.flatten_dict(v, new_key))
            else:
                items[new_key] = v
        return items

    def write_to_file(self, file_path):
            with open(file_path, 'w') as fd:
                input_data_str = []
                for key, value in self.input.items():
                    # Use f-strings for better readability
                    # Use json.dumps to format the value nicely if it's a dictionary or list
                    formatted_value = json.dumps(value, indent=2) if isinstance(value, (dict, list)) else str(value)
                    input_data_str.append(f"{key} {{ {formatted_value} }}\n")
                fd.write( ''.join(input_data_str ))

                # self._legacy_styling(fd)

    def _legacy_styling(self,fd):
            self.flat_input=self.matter.params
            input_data_str = []

            for key, value in self.flat_input.items():
                # Use f-strings for better readability
                # Use json.dumps to format the value nicely if it's a dictionary or list
                formatted_value = json.dumps(value, indent=2) if isinstance(value, (dict, list)) else str(value)
                input_data_str.append(f"{key} {{ {formatted_value} }}\n")
            fd.write( ''.join(input_data_str ))
            species = sorted(set(self.ocean_structure.numbers))
            fd.write('znucl {{ {} }}\n'.format(' '.join(str(Z) for Z in species)))
            fd.write('typat')
            fd.write('{\n')
            types = []
            for Z in self.atoms.numbers:
                for n, Zs in enumerate(species):
                    if Z == Zs:
                        types.append(n + 1)
            n_entries_int = 20  # integer entries per line
            for n, type in enumerate(types):
                fd.write(' %d' % (type))
                if n > 1 and ((n % n_entries_int) == 1):
                    fd.write('\n')
            fd.write(' }\n')

            atomic_positions_str = []
            for atom in self.atoms:
                atomic_positions_str.append( '{coords[0]:.10f} {coords[1]:.10f} {coords[2]:.10f}\n'.format(
                    coords=[atom.a, atom.b, atom.c] ))

            fd.write( 'xred {\n' )
            fd.write( ''.join(atomic_positions_str))
            fd.write( '}\n' )

            fd.write( 'acell {{ {acell[0]} {acell[0]} {acell[0]} }} \n'.format( acell=[1/Bohr] ) )

            fd.write( 'rprim {{ {cell[0][0]:.14f} {cell[0][1]:.14f} {cell[0][2]:.14f}\n'
                    '        {cell[1][0]:.14f} {cell[1][1]:.14f} {cell[1][2]:.14f}\n'
                    '        {cell[2][0]:.14f} {cell[2][1]:.14f} {cell[2][2]:.14f}  }}\n'
                        ''.format(cell=self.atoms.cell))
    
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