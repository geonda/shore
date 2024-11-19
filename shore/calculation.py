import os
import subprocess
import time
import shutil
import json
from pathlib import Path
import re


import numpy as np
from tqdm import tqdm, tqdm_notebook
import matplotlib.pyplot as plt

import ase.io
from ase.atoms import Atoms
from ase.visualize import view as avw
from ase.build import make_supercell  # Assuming you may need this function later
from ase.units import Bohr
from shore.remote import RemoteServerManager
import pickle
import spglib

from IPython.display import JSON

# Shore-related imports
# from shore.input_manager import input_processing
import shore.jobcreator as jc
from shore import gview
from shore.workflow_graph import workflow as gworkflow
from shore.input_manager import *
from shore.results import ResultsHandler
import copy
# Optional: If you plan to use threading, you can keep this import.
import threading






class Calculation():
    def __init__(self, server=None, workflow=None, input=None):
        if workflow:
            self.workflow=workflow.workflow
            if workflow.server:
                self.server=workflow.server
                self.root=workflow.root
        elif server:
            self.server=server
            self.workflow=gworkflow()
            self.root=os.getcwd()
        else:
            self.server=None
            self.workflow=gworkflow()
            self.root=os.getcwd()
        self.fresh=True
        os.chdir(self.root)


        
        self.structure=input.structure
        self.name=input.name
        self.element = input.content.element

        self.rpath=f'{self.structure.name}/{self.name}'
        self.edge = input.content.edge
        self.workflow = self.workflow
        self.input=input
        self.local_dir=f'{self.root}/{self.structure.name}/{self.name}'
        self.log_dir=f'{self.root}/{self.structure.name}/{self.name}/logs'
        self.res_dir=f'{self.local_dir}/results/'
        self.path_to_jar=os.path.join(self.root,'jar')
        self.filename=f'{self.log_dir}/log_{input.name}_state'
        self.input_data=input.content.input
        self.job_id=None
        self.res=None



        self.workflow.add_node(str(self.structure.name), layer='structure', state='active')
        self.workflow.add_instance(node1=self.structure.name, node2=self.name, layer='input', )
        self.workflow.graph.nodes(data=True)

        self.list_ids=self._get_element_order(self.element)

        Path(self.local_dir).mkdir(parents=True, exist_ok=True)
        Path(self.log_dir).mkdir(parents=True, exist_ok=True)
        Path(self.res_dir).mkdir(parents=True, exist_ok=True)

        # create ocean.in
        self.input.content.write_to_file(f'{self.local_dir}/ocean.in')
        self.input.light.write_to_folder(file_path=f'{self.local_dir}')
        
        if self.server:
            self.remote_dir=f'{self.server.root}/{self.rpath}'
            self.server.connect()
            self.server.remote_dir_init(f"{self.rpath}")
            self.server.upload_file(f"{self.local_dir}/ocean.in", self.remote_dir)
            for id,_ in enumerate(self.input.light.photons):
                self.server.connect()
                self.server.upload_file(f"{self.local_dir}/photon{id+1}", self.remote_dir)
            jc.JobScriptCreator(ncores=self.server.cores).generate_script(path=self.local_dir, command=self.server.command)
            self.server.connect()
            self.server.upload_file(f"{self.local_dir}/job.sh", self.remote_dir)
     

        self.stages=dict(
                        parsing=['Storing parsed data','Finished running extractPsp','Done with parsing'],
                        atomic=['Entering OPF stage','Entering DFT stage'],
                        dft=['Entering DFT stage','DFT for BSE final states complete','DFT section is complete'],
                        prep=['Entering PREP stage','Entering SCREENing stage'], 
                        screen=['Entering SCREENing stage','Time offset:'],
                        bse=['CNBSE stage','Ocean is done']
                        )
        
        self.stages_states=dict(parsing=0, opf=0, dft=0, prep=0, screen=0, cnbse=0)
        self.workflow.graph.nodes[f'{self.name}']['state']='active'
        self.workflow.add_instance(node1=f'{self.name}',node2=f'{self.name}-parsing', layer='parsing', )
        self.workflow.add_instance(node1=f'{self.name}-parsing',node2=f'{self.name}-atomic', layer='atomic', )
        self.workflow.add_instance(node1=f'{self.name}-atomic',node2=f'{self.name}-dft', layer='dft', )   
        self.workflow.add_instance(node1=f'{self.name}-dft',node2=f'{self.name}-prep', layer='prep', )   
        self.workflow.add_instance(node1=f'{self.name}-prep',node2=f'{self.name}-screen', layer='screen', )   
        self.workflow.add_instance(node1=f'{self.name}-screen',node2=f'{self.name}-bse', layer='bse', ) 
        self.workflow.add_instance(node1=f'{self.name}-bse', 
                                         node2=f'{self.name}-results', layer='results')
        for atom_sites in self.list_ids:
            self.workflow.add_instance(node1=f'{self.name}-results',
                                            node2=f'{self.name}-{self.edge}-{self.element}-{atom_sites}',
                                            layer='xas results',)
        
        self._check_if_there_is_something()
        if not self.fresh:
            self.load_attributes()

  
        
    def _get_element_order(self,element):
        id=1
        output=[]
        for s,n  in zip(self.structure.atoms.get_chemical_symbols(),self.structure.atoms.get_atomic_numbers()):
            if s==element:
                output.append(id)
                id+=1
        return output
    
    def set(self, key, value):
        """Update the value of a specific key in the dictionary and save if changed."""
        if key in self.input_data:
            if self.input_data[key] != value:  # Check if the value has changed
                print(f'Changing {key} from {self.input[key]} to {value}')
                self.input_data[key] = value
                self.save_attributes()  # Save when changed
        else:
            print(f'Adding new attribute {key} with value {value}')
            self.input_data[key] = value
            self.save_attributes() 

    def load_attributes(self):
        """Load attribute values from a JSON file."""
        if os.path.exists(self.filename):
            with open(self.filename, 'r') as fin:
                attributes = json.load(fin)
                for key, value in attributes.items():
                    setattr(self, key, value)  # Set each attribute dynamically


    def is_serializable(self, obj):
        """Check if an object is serializable."""
        return isinstance(obj, (int, float, str, list, dict))

    def save_attributes(self):
        """Save current attribute values to a JSON file."""
        attributes = vars(self)  # Automatically get all attributes
        
        # Filter out non-serializable attributes
        serializable_attributes = {
            key: value for key, value in attributes.items() if self.is_serializable(value)
        }
        
        with open(self.filename, 'w') as fout:
            json.dump(serializable_attributes, fout, indent=4)
    
    def handle_input(self):
        self.input.content.write_to_file(f"{self.local_dir}/ocean.in")
        
        if self.server:
            self.server.connect()
            self.server.upload_file(f"{self.local_dir}/ocean.in", self.remote_dir)
            
    def _check_if_there_is_something(self):
        
        for k,v in self.stages_states.items():
            if v>0:
                self.workflow.graph.nodes[f'{self.name}-{k}']['state']='active'
            # else:
            #     print('No')
    
    def sync(self):
        # self.get_state()

        # Prepare to download files and count total files
        total_files = 0

        # Count files for spectra download
        try:
            spectra_files = self.server.list_files(f'{self.remote_dir}/CNBSE/')  # Assuming this method lists files
            total_files += len(spectra_files)
        except Exception as e:
            print(e)

        # Count files for DFT download
        dft_files = ['scf.in', 'scf.out']
        total_files += len(dft_files)

        # Count log files
        log_files = ['mpi_avg.log', 'ocean.log']
        total_files += len(log_files)

        # Create a single progress bar for all downloads
        with tqdm(total=total_files, desc='Syncing Files', unit='file') as pbar:
            # Download spectra files
            try:
                for file_name in spectra_files:
                    self.server.download_file(file_name, f'{self.local_dir}/results/', f'{self.remote_dir}/CNBSE/')
                    pbar.update(1)  # Update progress bar after each file download
            except Exception as e:
                print(e)

            # Download DFT files
            try:
                for name in dft_files:
                    self.server.download_file(name, f'{self.local_dir}/logs', f'{self.remote_dir}/DFT/')
                    pbar.update(1)
            except Exception as e:
                print(e)

            # Download log files
            try:
                Path(f'{self.local_dir}/logs').mkdir(parents=True, exist_ok=True)
                for name in log_files:
                    if name == 'mpi_avg.log':
                        self.server.download_file(name, f'{self.local_dir}/logs', f'{self.remote_dir}/SCREEN/')
                    elif name == 'ocean.log':
                        self.server.download_file(name, f'{self.local_dir}/logs', f'{self.remote_dir}/CNBSE/')
                    pbar.update(1)
            except Exception as e:
                print(e)

        self.stages_states['results']=1
        self.workflow.graph.nodes[f'{self.name}-results']['state']='active'
        for atom_sites in self.list_ids:
            self.stages_states[f'{self.name}-{self.edge}-{self.element}-{atom_sites}']=1
        self.workflow.graph.nodes[f'{self.name}-{self.edge}-{self.element}-{atom_sites}']['state']='active'

        self.res = ResultsHandler(path=f'{self.local_dir}/results/', name=self.name)
        return self.res
                
    def __sync(self,):
            self.get_state()
            try:
                self.server.download_spectra( f'{self.remote_dir}/CNBSE/',f'{self.local_dir}/results/',)
            except Exception as e:
                print(e)
            try:
                for name in ['scf.in', 'scf.out', 'nscf.in', 'nscf.out']:
                    self.server.download_file(name, f'{self.local_dir}/dft',f'{self.remote_dir}/DFT/',)
            except Exception as e:
                print(e)
            # try:
            #     self.server.download_folder(f'{self.remote_dir}/Common/',f'{self.local_dir}/logs/common',)
            # except Exception as e:
            #     print(e)
            try:
                Path(f'{self.local_dir}/logs').mkdir(parents=True, exist_ok=True)
                for name in ['mpi_avg.log']:
                    self.server.download_file(name, f'{self.local_dir}/logs',f'{self.remote_dir}/SCREEN/',)
            except Exception as e:
                print(e)

            try:
                Path(f'{self.local_dir}/logs').mkdir(parents=True, exist_ok=True)
                for name in ['ocean.log']:
                    self.server.download_file(name, f'{self.local_dir}/logs',f'{self.remote_dir}/CNBSE/',)
            except Exception as e:
                print(e)
            self.res=ResultsHandler(path=f'{self.local_dir}/results/', name=self.name)
            return self.res
      

    def _edge_short(self, edge):
        if edge=='K':
            return '1s'
        else:
            return '2p'
        

    def _check_folder_exists_with_content(self,folder_path):
        """
        Check if a folder exists and contains files.

        Parameters:
        folder_path (str): The path to the folder to check.

        Returns:
        bool: True if the folder exists and contains files, False otherwise.
        """
        # Check if the folder exists
        if os.path.isdir(folder_path):
            # Check if the folder contains any files
            if any(os.path.isfile(os.path.join(folder_path, f)) or os.path.isdir(os.path.join(folder_path, f)) for f in os.listdir(folder_path)):
                return True  # Folder exists and has files
            else:
                return False  # Folder exists but is empty
        else:
            return False  # Folder does not exist
        

    def run(self, overwrite=False,monitor=False):
        if self.server:
            self.save_attributes()
            self._run_remote(overwrite=overwrite, monitor=monitor)
        else:
            self._run_local(overwrite=overwrite, monitor=monitor)
        

    
    def _run_local(self, overwrite=False, monitor=True):
        """Execute a bash command in the directory under path."""

        if self.path is None:
            raise ValueError("Path is not set. Please provide a valid path.")
        
        # Change the current working directory to the specified path
        try:
            os.chdir(self.path)  # Change to the specified directory
            self.handle_input()


            if self._check_folder_exists_with_content(f'{self.local_dir}/CNBSE/') and overwrite==False:
                print('Heavy part is already done. if you want to rerun it use overwrite=true')
                # for atom_sites in self.list_ids:
                #     self.workflow.add_instance(node1=f'{self.element}-{self.edge} edge',
                #                             node2=f'XAS {self.element} {atom_sites}', layer='xas results',)
            else:                          
                os.system('rm -r ./*')
                os.system(f'cp {self.root}/ocean.sh {self.local_dir}/')
                # print(f"Changed directory to: {self.path}")
                #
                # result = subprocess.run(["./ocean.sh"], capture_output=True, text=True)
                result = subprocess.Popen('./ocean.sh', shell=True, text=True, 
                                          stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE)
                if monitor:
                    self.monitor()

        except FileNotFoundError:
            print(f"Error: The directory {self.local_dir} does not exist.")
        
    def _run_remote(self, overwrite=False, monitor=True):
        
        try:
            self.handle_input()
            self.server.connect()
            if self.server.check_folder_exists_and_not_empty(f'{self.remote_dir}/CNBSE/') and overwrite==False:
                print('Heavy part is already done. if you want to rerun it use overwrite=true')
            else:  
                if self.server.sbatch:                      
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
                else:

                    # command=f'source ~/miniforge3/bin/activate new ; cd {self.remote_dir}; /home/a.geondzhian/bin/ocean-acbn0/ocean.pl ocean.in > log &'
                    command=f'source /etc/profile.d/modules.sh ; module load q-ch/qe/7.3.1/gcc/11.2/mpich/mkl; cd {self.remote_dir}; /home/a.geondzhian/bin/ocean-acbn0/ocean.pl ocean.in > log &'
                    # command=f'cd {self.remote_dir}; pwd; /home/a.geondzhian/bin/ocean-acbn0/ocean.pl ocean.in > log'
                    transport=self.server.ssh_client.get_transport()
                    channel=transport.open_session()
                    try:
                        # Execute the command
                        channel.exec_command(command)
                    finally:
                        # Close the channel to free resources
                        channel.close()
                    
                
                if monitor:
                    self._remote_monitor()

        except Exception as e:
            print(f"{e}")

    def __read_error_file(self,file_path):
        """
        Reads the contents of an error file into a variable.

        :param file_path: The path to the error file.
        :return: The contents of the error file as a string.
        """
        try:
            with open(file_path, 'r') as file:
                error_contents = file.read()
            return error_contents
        except FileNotFoundError:
            print(f"Error: The file '{file_path}' was not found.")
            return None
        except Exception as e:
            print(f"An error occurred while reading the file: {e}")
            return None
    
    def cancel(self, job_id=None):
        if job_id:
            tmp_id=job_id
        elif self.job_id:
            tmp_id=self.job_id
        else:
            print("no id found")
        try:
            self.server.connect()
            stdin, stdout, stderr=self.server.ssh_client.exec_command(f'scancel {tmp_id}')
            output = stdout.read().decode('utf-8')
            error_output = stderr.read().decode('utf-8')
            print(output,error_output)
        except:
            print("didn't go through")
            pass
            
        
    def get_status(self):
        self.get_state()
        print('Errors:')
        for line in self.err.split('\n'):
            print(line)
        print('Messages:')
        for line in self.out.split('\n'):
            print(line)
        print('Now is runing:')
        try:
            self.server.connect()
            stdin, stdout, stderr=self.server.ssh_client.exec_command(f'squeue --user={self.server.username}')
            output = stdout.read().decode('utf-8')
            error_output = stderr.read().decode('utf-8')
            print(output,error_output)
        except:
            print('Nothing')


    def get_state(self):
        self.file_path=f'{self.local_dir}/logs/log'
        self.server.connect()
        self.server.download_file('log',f'{self.local_dir}/logs',self.remote_dir)
        server_updates=False
        try:
            self.server.download_file('err',os.path.join(self.local_dir,'logs'),self.remote_dir)
            self.server.download_file('out',os.path.join(self.local_dir,'logs'),self.remote_dir)
            server_updates=True
        except:
            print('no info from the slurm')
        self.server.disconnect()
        if server_updates:
            self.err=self.__read_error_file(f'{self.local_dir}/logs/err')
            self.out=self.__read_error_file(f'{self.local_dir}/logs/out')
    
        

        for stage,desc in self.stages.items():
            count=0
            with tqdm_notebook(total=len(desc),desc=f'{stage}') as pbar:
                for steps in desc:
                    if  self._progress_info(f'{self.local_dir}/logs/log',steps):
                        pbar.update(1) 
                        count+=1
            if count ==len(desc):
                self.stages_states[stage]=1
                self.workflow.graph.nodes[f'{self.name}-{stage}']['state']='active'


        self.save_attributes()
        # print("Process completed.")


    def _remote_monitor(self):
        self.progress = {key: 0 for key in ['initial', 'opf', 'dft', 'prep', 'screen', 'cnbse', 'ocean']}
        self.file_path = f'{self.local_dir}/logs/log'
        
        # Flag to control the monitoring thread
        self.monitoring_active = True

        try:
            self.server.connect()
            self.monitor_thread = threading.Thread(target=self.server.monitor_files,
                                                        args=(self.remote_dir, ['log'],
                                                        f"{self.local_dir}/logs",self.monitoring_active))
            self.monitor_thread.start()

            with tqdm(total=100, desc='Ocean is rising') as pbar:
                while self.progress['ocean'] < 100 and self.monitoring_active:
                    time.sleep(1)
                    self.extract_progress_info(self.file_path)
                    max_progress = max(self.progress.values())
                    pbar.update(max(0, max_progress - pbar.n))

        except KeyboardInterrupt:
            print("Process was interrupted by the user.")
            self.monitoring_active = False  # Signal the thread to stop
        except Exception as e:
            print(f"An error occurred: {e}")
            self.monitoring_active = False  # Signal the thread to stop
        finally:
            # print("Cleaning up...")
            self.server.disconnect()

        # Update workflow state and save attributes
        self.workflow.graph.nodes[f'{self.name}']['state'] = 'active'
        self.save_attributes()
        # print("Process completed.")



    def monitor_remote(self):

        for folder,files in self.dict_to_check.items():
            if not os.path.isdir(f"{self.local_dir}/{folder}"):
                os.mkdir(f"{self.local_dir}/{folder}")
            
            self.server.connect()
            self.monitor_thread = threading.Thread(target=self.server.monitor_files,
                                                    args=(f"{self.remote_dir}/{folder}", files,
                                                            f"{self.local_dir}/{folder}"))
            self.monitor_thread.start()
            time.sleep(10)
            with tqdm_notebook(total=3,desc=f'{folder} stage is running') as pbar:
                while not self._check_files_exist(f"{self.local_dir}/{folder}",files):
                    time.sleep(1) 
                    for file in files:
                        if os.path.isfile(os.path.join(f"{self.local_dir}/{folder}", file)):
                            pbar.update(1) 
            # self.monitor_thread.join()
            # self.server.disconnect()

            self.workflow.graph.nodes[f'{self.name}']['state']='active'
        # self.workflow.graph.nodes[f'{self.element}']['state']='active'
        # self.workflow.graph.nodes[f'{self.element}-{self.edge} edge']['state']='active'
        self.save_attributes()
        print("Process completed.")

    def _progress_info(self, file, find):
        with open(file, 'r') as file:
            for line in file:
                if find in line:
                    return True
                    
        return False
    
    def extract_progress_info(self, file):
        with open(file, 'r') as file:
            for line in file:
                if 'Welcome to OCEAN' in line:
                    self.progress['initial'] = 1
                if 'Entering OPF stage' in line:
                    self.progress['opf'] = 10
                if 'Entering DFT stage' in line:
                    self.progress['dft'] = 30
                if 'SCF stage complete' in line:
                    self.progress['scf'] = 40
                if 'Entering PREP stage' in  line:
                    self.progress['prep'] = 50
                if 'Entering SCREENing stage' in line:
                    self.progress['screen']=60
                if 'Entering CNBSE stage' in line:
                    self.progress['cnbse']=80
                if 'Ocean is done' in line: 
                    self.progress['ocean']=100

    def monitor(self):
        self.progress = {
            'initial': 0,
            'opf': 0,
            'dft': 0,
            'prep': 0,
            'screen': 0,
            'cnbse': 0,
            'ocean': 0,
        }
        time.sleep(5)
        self.file_path='/home/geonda/shore/log'

        with tqdm_notebook(total=100,desc=f'ocean is rising') as pbar:
            try:
                while not self.progress['ocean']==100:
                    time.sleep(1) 
                    self.extract_progress_info()
                    max_progress=max(self.progress.values())
                    pbar.update(max_progress-pbar.n) 
                        
            except KeyboardInterrupt:
                    # Handle the case where the user interrupts the process
                    print("Process was interrupted by the user.")
            finally:
                pbar.close()  # Ensure the progress bar is closed properly
                # print("Cleaning up...")
                
        self.workflow.graph.nodes[f'{self.name}']['state']='active'
        for item,value in self.dict_to_check.items():
            if value:
                self.workflow.graph.nodes[f'{self.element}_{self.edge}_{item}']['state']='active'
        self.save_attributes()
        print("Process completed.")
    def _check_files_exist(self,folder,files):
        return all(os.path.isfile(os.path.join(folder, file)) for file in files)
    def _adjust_zeros(self, number):
        if number > 9:
            return str(number)
        else:
            return f'0{number}'
    def _adjust_element(self,element):
        if self.edge=='L':
            return element
        else:
            return f"{element}"
        
    def plot(self,sites=[],pols=[], label=None, energy_shift=None):
        self.x=0
        self.y=0
        _sites= sites if sites!=[] else self.list_ids
        _pols = pols if pols!=[] else [1,2,3]
        for site in _sites:
            for pol in _pols:
                try:
                    data=np.loadtxt(f'{self.local_dir}/results//absspct_{self._adjust_element(self.element)}.00{self._adjust_zeros(site)}_{self._edge_short(self.edge)}_0{pol}').transpose()
                    self.x=data[0]
                    self.y+=data[1]
                    self.workflow.graph.nodes[f'{self.name}-{self.edge}-{self.element}-{site}']['state']='active'
                except Exception as e:
                    print(f'Probably no file {e}')
                    self.workflow.graph.nodes[f'{self.name}-{self.edge}-{self.element}-{site}']['state']='inactive'
        try:
            if energy_shift:
                plt.plot(self.x+energy_shift,self.y, label=label)
            else:
                plt.plot(self.x,self.y, label=label)
            self.save_attributes()
        except Exception as e:
            print(f'Probably no file {e}')
        plt.ylabel('XAS Intensity, arb. units')
        plt.xlabel('Energy, eV')
    
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
    
    def clean_local_dir(self):
        """Delete the local directory."""
        if os.path.exists(self.local_dir):
            shutil.rmtree(self.local_dir)
            print(f"Local directory {self.local_dir} has been deleted.")
        else:
            print(f"Local directory {self.local_dir} does not exist.")

    def clean_remote_dir(self):
        """Delete the remote directory."""
        if self.server:
            self.server.connect()
            try:
                self.server.delete_directory(self.remote_dir)  # Assuming delete_directory is a method in your server class
                print(f"Remote directory {self.remote_dir} has been deleted.")
            except Exception as e:
                print(f"Error deleting remote directory: {e}")
            finally:
                self.server.disconnect()
        else:
            print("No server connection available to delete the remote directory.")
    
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