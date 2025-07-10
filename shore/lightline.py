import fdmnes 
import numpy as np
import time
import sys
import matplotlib.pyplot as plt
import os
import glob
import numpy as np
from collections import defaultdict
import shore.jobcreator as jc
import re

class XASSpectraRunner:
    def __init__(self,id=0,dir='./', server=None, cif_path=None,  resonant_atom=None, fdmnes_path='/opt/api/fdmnes_linux64'):
        self.local_dir=dir
        os.makedirs(self.local_dir, exist_ok=True)
        self.cif_path = cif_path
        self.resonant_atom = resonant_atom
        self.fdmnes_path=fdmnes_path
        self.current_id=id
        self.sim=fdmnes.fdmnes(self.cif_path, resonant=self.resonant_atom, 
                          fdmnes_path=self.fdmnes_path)
        self.sim.P.Range = (-5, 0.1,5,0.5, 20)
        self.sim.P.radius = 6.0
        self.sim.P.Quadrupole = True
        self.sim.P.Convolution = True
        self.sim.P.Rpotmax = 13.
        self.sim.P.Green = True

        self.server=server
        self.input_file=f'xas_inp_{self.current_id}'
        self.input_name=f'{self.local_dir}/xas_inp_{self.current_id}'
        self.sim.WriteInputFile(self.input_name, overwrite=True)
        self.init_fdmnes()
        
        self.rpath=f'matsolver_workdir/{self.current_id}'
        self.remote_dir=f'{self.server.root}/{self.rpath}'

    def init_fdmnes(self):
        output = []
        output.append("1")
        output.append("")
        output.append(self.input_file)
        with open(f"{self.local_dir}/fdmfile.txt", "wb") as f:
            f.write(os.linesep.join(output).encode())

    def run(self):
            
        self.server.connect()
        self.server.remote_dir_init(f"{self.rpath}")
        self.server.upload_file(self.input_name, self.remote_dir)
        jc.JobScriptCreator(ncores=self.server.cores).generate_script(path=self.local_dir, command=self.server.command)
        self.server.connect()
        self.server.upload_file(f"{self.local_dir}/job.sh", self.remote_dir)
        self.server.connect()
        self.server.upload_file(f"{self.local_dir}/fdmfile.txt", self.remote_dir)
        self._run_remote()
    
    def check_status(self):
        self.status='Running'
        try:
            self.server.download_file('out', f'{self.local_dir}/', f'{self.remote_dir}/') 
            try: 
                with open(f'{self.local_dir}/out') as f:
                    lines=f.readlines()
                    for line in lines:
                        if 'Arctangent' in line:
                            self.status='FINISHED'
            except:
                self.status='Running'
        except:
            print('Can not access output file')
    
    def get_spectra(self):
        try:
            self.server.download_file(f'xas_out_{self.current_id}.txt', f'{self.local_dir}/', f'{self.remote_dir}/') 
        except:
            print('error reading xas')
        

        
    
    def _run_remote(self,):
        
        self.server.connect()
                        
        command=f'cd {self.remote_dir}; pwd; sbatch job.sh'
        stdin, stdout, stderr=self.server.ssh_client.exec_command(command)
                    # transport.close()
        output = stdout.read().decode('utf-8')
        error_output = stderr.read().decode('utf-8')
        if error_output:
            print(f"Error submitting job: {error_output}")
        match = re.search(r'Submitted batch job (\d+)', output)
        if match:
            job_id = match.group(1)
            print(f"Job submitted successfully with Job ID: {job_id}")
            self.job_id=job_id
        else:
            print("Could not retrieve Job ID from sbatch output.")


    def read(self, conv=True):
        # Pattern for numbered files, e.g. xas_out_123_1.txt, xas_out_123_2.txt ...
        if not conv:
            numbered_pattern = os.path.join(self.local_dir, f'xas_out_{self.current_id}_*.txt')
            all_files = sorted(glob.glob(numbered_pattern))

            # Exclude files ending with _bav.txt or _conv.txt
            numbered_files = [
                f for f in all_files
                if not (f.endswith('_bav.txt') or f.endswith('_conv.txt'))
            ]

            count = 0
            spectra = {}

            if numbered_files:
                # Multiple atom spectra files found            
                for file in numbered_files:
                    # Load data skipping first 2 rows (adjust if needed)
                    data = np.loadtxt(file, skiprows=2).transpose()
                    energy=self.read_first_element_of_first_line(file)
                    data[0]+=energy
                    spectra[count] = data.tolist()
                    count += 1
                return spectra  # Dict of spectra, one per atom
            else:
                # Fallback to single file
                single_file = os.path.join(self.local_dir, f'xas_out_{self.current_id}.txt')
                if os.path.exists(single_file):
                    data = np.loadtxt(single_file, skiprows=2).transpose()
                    energy=self.read_first_element_of_first_line(file)
                    data[0]+=energy
                    spectra[count] = data.tolist()
                    return spectra
                else:
                    raise FileNotFoundError(f"No XAS spectrum files found for id {self.current_id}")
        else:
            count = 0
            spectra = {}
            single_file = os.path.join(self.local_dir, f'xas_out_{self.current_id}_conv.txt')
            if os.path.exists(single_file):
                data = np.loadtxt(single_file, skiprows=2).transpose()
                # energy=self.read_first_element_of_first_line(single_file)
                energy=529.
                data[0]+=energy
                spectra[count] = data.tolist()
                return spectra
            
    def read_first_element_of_first_line(self,filename):
        with open(filename, 'r') as f:
            first_line = f.readline().strip()  # read first line and strip whitespace
            first_element = first_line.split()[0]  # split by whitespace and take first element
            return float(first_element)

# if __name__=='__main__':
#     test=XASSpectraRunner(id=0,dir='./test', cif_path='./BaTiO3.cif', resonant_atom='O', )
#     # test.run()
#     test.read()