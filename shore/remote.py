import paramiko
import time
import threading
import logging
import os
import pickle
from pathlib import Path

class RemoteServerManager:
    def __init__(self, remote_host=None, username=None,key=None, password=None, load=None, command='/opt/cms/q-ch/ocean/ocean.pl ocean.in > log'):
        if not load:
            self.remote_host = remote_host
            self.username = username
            self.password = password
            if not key:
                self.file_name=f'/home/{username}/.ssh/id_rsa'
            else:
                self.file_name=key
            # self.root=f'/home/{username}/ocean_work_dir/'
            self.root=self._init_directory_if_not_exists('/ocean_work_dir')
            self.sbatch=True
            self.monitor_active=False
            self.command=command
        else:
            self.load(load)

        self.cores=4
        self.server_dir='./server_logs/'
        Path(self.server_dir).mkdir(parents=True, exist_ok=True)
        # Set up logging
        logging.basicConfig(
            filename=f'{self.server_dir}/remote_server_manager.log',  # Log file name
            level=logging.INFO,                     # Log level
            format='%(asctime)s - %(levelname)s - %(message)s'  # Log format
        )
        
        self.ssh_client = None
        self.sftp_client = None
        self.remote_dir=None

    def save(self, path):
        with open(f'{path}', 'wb') as f:
            pickle.dump(self.__dict__,f)
    
    def load(self,path):
        with open(f'{path}', 'rb') as f:
            tmp=pickle.load(f)
        self.__dict__.update(tmp)

    def delete_directory(self, remote_dir):
        """Delete a directory on the remote server.

        Args:
            remote_dir (str): The path to the remote directory to be deleted.
        """
        command = f'rm -rf {remote_dir}'  # Command to remove directory and its contents
        stdin, stdout, stderr = self.ssh_client.exec_command(command)

        # Capture output and errors
        output = stdout.read().decode('utf-8')
        error_output = stderr.read().decode('utf-8')

        if error_output:
            print(f"Error deleting directory: {error_output}")
        else:
            print(f"Successfully deleted directory: {remote_dir}")
        
        
    

    def remote_dir_init(self,rpath):
        self.create_remote_directory_if_not_exists(f"{rpath}")
        self.remote_dir=f"{self.root}/{rpath}"
    
    def _init_directory_if_not_exists(self,remote_path):
        """Recursively create a directory on remote server if not exists."""
        # Split the path into parts
        self.connect()
        parts = remote_path.split('/')
        
        path_to_create = f"/home/{self.username}/"
        for part in parts:
            if part:  # Ignore empty parts (in case of leading slash)
                path_to_create = os.path.join(path_to_create, part)
                try:
                    remote_dir_check = self.ssh_client.exec_command(f"if [ ! -d '{path_to_create}' ]; then echo 'not_exists'; fi")[1].read().strip().decode()
                    if remote_dir_check == 'not_exists':
                        self.ssh_client.exec_command(f"mkdir '{path_to_create}'")
                        logging.info(f"Created directory: {path_to_create}")
                    else:
                        logging.info(f"Directory already exists: {path_to_create}")
                except Exception as e:
                    logging.error(f"Failed to create directory '{path_to_create}': {e}")
        # self.disconnect()
        return f"/home/{self.username}/{remote_path}"
        
    
    def create_remote_directory_if_not_exists(self,remote_path):
        """Recursively create a directory on remote server if not exists."""
        # Split the path into parts
        parts = remote_path.split('/')
        
        path_to_create = self.root
        for part in parts:
            if part:  # Ignore empty parts (in case of leading slash)
                path_to_create = os.path.join(path_to_create, part)
                try:
                    remote_dir_check = self.ssh_client.exec_command(f"if [ ! -d '{path_to_create}' ]; then echo 'not_exists'; fi")[1].read().strip().decode()
                    if remote_dir_check == 'not_exists':
                        self.ssh_client.exec_command(f"mkdir '{path_to_create}'")
                        logging.info(f"Created directory: {path_to_create}")
                    else:
                        logging.info(f"Directory already exists: {path_to_create}")
                except Exception as e:
                    logging.error(f"Failed to create directory '{path_to_create}': {e}")



    def download_folder(self, remote_directory, local_directory):
        """
        Downloads a folder from a remote server to a local directory using Paramiko.

        :param ssh_client: An active SSH client connection.
        :param remote_directory: The path to the folder on the remote server.
        :param local_directory: The path to the local directory where files should be saved.
        """
        self.connect()
        sftp = self.ssh_client.open_sftp()
        
        try:
            # Create local directory if it doesn't exist
            if not os.path.exists(local_directory):
                os.makedirs(local_directory)
            
            # List files in the remote directory
            for item in sftp.listdir(remote_directory):
                remote_path = os.path.join(remote_directory, item)
                local_path = os.path.join(local_directory, item)
                
                # Check if the remote path is a directory
               
                sftp.get(remote_path, local_path)
                # print(f"Downloaded: {local_path}")
        finally:
            sftp.close()
    def list_files(self, remote_directory):
        """
        Downloads a folder from a remote server to a local directory using Paramiko.

        :param ssh_client: An active SSH client connection.
        :param remote_directory: The path to the folder on the remote server.
        :param local_directory: The path to the local directory where files should be saved.
        """
        self.connect()
        sftp = self.ssh_client.open_sftp()
        
        try:
         
            
            remote_files = sftp.listdir(remote_directory)

            # Filter files that start with 'absspct'
            absspct_files = [f for f in remote_files if f.startswith("absspct")]
            
            # List files in the remote directory
        except: 
            print('No spectra found')
            return None
        finally:
            sftp.close()
            return absspct_files
    
    def download_spectra(self, remote_directory, local_directory):
        """
        Downloads a folder from a remote server to a local directory using Paramiko.

        :param ssh_client: An active SSH client connection.
        :param remote_directory: The path to the folder on the remote server.
        :param local_directory: The path to the local directory where files should be saved.
        """
        self.connect()
        sftp = self.ssh_client.open_sftp()
        
        try:
            # Create local directory if it doesn't exist
            if not os.path.exists(local_directory):
                os.makedirs(local_directory)
            
            remote_files = sftp.listdir(remote_directory)

            # Filter files that start with 'absspct'
            absspct_files = [f for f in remote_files if f.startswith("absspct")]
            print(absspct_files)
            # List files in the remote directory
            for item in absspct_files:
                remote_path = os.path.join(remote_directory, item)
                local_path = os.path.join(local_directory, item)
                
                # Check if the remote path is a directory
               
                sftp.get(remote_path, local_path)
                # print(f"Downloaded: {local_path}")
        finally:
            sftp.close()

    # def is_remote_directory(self, sftp, remote_path):
    #     """
    #     Check if the remote path is a directory.

    #     :param sftp: An active SFTP session.
    #     :param remote_path: The path to check.
    #     :return: True if the path is a directory, False otherwise.
    #     """
    #     try:
    #         return stat.S_ISDIR(sftp.stat(remote_path).st_mode)
    #     except IOError:
    #         return False
    

    def connect(self):
        """Establish an SSH connection to the remote server."""
        self.ssh_client = paramiko.SSHClient()
        self.ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
      
        try:
            private_key = paramiko.RSAKey.from_private_key_file(self.file_name)
            self.ssh_client.connect(self.remote_host, username=self.username, password=self.password, pkey=private_key)
            self.sftp_client = self.ssh_client.open_sftp()
            logging.info(f"Connected to {self.remote_host}")
        except Exception as e:
            logging.error(f"Failed to connect to {self.remote_host}: {e}")
            raise

    def disconnect(self):
        """Close the SSH connection."""
        if self.sftp_client:
            self.sftp_client.close()
            logging.info(f"SFTP session closed.")
        if self.ssh_client:
            self.ssh_client.close()
            logging.info(f"Disconnected from {self.remote_host}")

    def create_folder(self, folder_path):
        """Create a new folder on the remote server."""
        command = f'mkdir -p {folder_path}'  # Use -p to avoid error if the folder already exists
        try:
            stdin, stdout, stderr = self.ssh_client.exec_command(command)
            exit_status = stdout.channel.recv_exit_status()  # Wait for command to finish
            
            if exit_status == 0:
                logging.info(f"Successfully created folder: {folder_path}")
                return True
            else:
                error_message = stderr.read().decode().strip()
                logging.error(f"Error creating folder {folder_path}: {error_message}")
                return False
        except Exception as e:
            logging.error(f"Exception occurred while creating folder: {e}")
            return False

    def upload_file(self, local_file_path, remote_directory):
        """Upload a file from local to the remote directory."""
        try:
            remote_file_path = f"{remote_directory}/{local_file_path.split('/')[-1]}"  # Get the filename from local path
            
            with self.sftp_client as sftp:
                sftp.put(local_file_path, remote_file_path)  # Upload the file
                logging.info(f"Uploaded {local_file_path} to {remote_file_path}")
                # print(f"Uploaded {local_file_path} to {remote_file_path}")
                
        except Exception as e:
            logging.error(f"Error uploading file {local_file_path}: {e}")

    def execute_command(self, command):
        """Execute a command on the remote server and return the output."""
        try:
            stdin, stdout, stderr = self.ssh_client.exec_command(command)
            exit_status = stdout.channel.recv_exit_status()  # Wait for command to finish
            
            output = stdout.read().decode().strip()
            error_output = stderr.read().decode().strip()

            if exit_status == 0:
                logging.info(f"Command executed successfully: {command}")
                return output
            else:
                logging.error(f"Error executing command '{command}': {error_output}")
                return None, error_output
        except Exception as e:
            logging.error(f"Exception occurred while executing command: {e}")
            return None
 
    
    def reconnect(self):
        try:
            self.disconnect()  # Ensure previous connection is closed
            self.connect()      # Reconnect to the server
            print("Reconnected successfully.")
        except Exception as e:
            print(f"Failed to reconnect: {e}")


    def monitor_files(self, remote_directory, files_to_monitor, local_directory, monitoring_active,polling_interval=5,):
        """Monitor specified files in the remote directory and download them if modified."""
        previous_state = {}
        
        while monitoring_active:
            current_state = {}
            
            for filename in files_to_monitor:
                try:
                    mtime = self.get_file_modification_time(filename, remote_directory)
                    if mtime is not None:
                        current_state[filename] = mtime
                        
                        # Check for new or modified files
                        if filename not in previous_state:
                            logging.info(f"New file detected: {filename}")
                            self.download_file(filename, local_directory, remote_directory)
                        elif mtime != previous_state[filename]:
                            logging.info(f"File modified: {filename}")
                            self.download_file(filename, local_directory, remote_directory)  # Download modified file
                    else:
                        if filename in previous_state:
                            logging.info(f"File deleted: {filename}")

                except OSError as e:
                    if "Socket is closed" in str(e):
                        logging.warning("Socket closed unexpectedly. Attempting to reconnect...")
                        self.reconnect()  # Implement a reconnect method
                    else:
                        logging.error(f"An OSError occurred while monitoring {filename}: {e}")
                except Exception as e:
                    logging.error(f"An unexpected error occurred while monitoring {filename}: {e}")

            previous_state = current_state
            time.sleep(polling_interval)  # Wait before polling again

   

    def get_file_modification_time(self, filename, remote_directory):
        """Retrieve the modification time of a specific file."""
        filepath = f"{remote_directory}/{filename}"
        try:
            file_info = self.sftp_client.stat(filepath)
            return file_info.st_mtime  # Return modification time
        except FileNotFoundError:
            return None  # File does not exist

    def download_file(self, filename, local_directory, remote_directory):
        """Download the specified file from the remote directory to the local directory."""
        remote_filepath = f"{remote_directory}/{filename}"
        local_filepath = f"{local_directory}/{filename}"
        
        try:
            logging.info(f"Downloading {filename}...")
            self.sftp_client.get(remote_filepath, local_filepath)  # Download the file
            logging.info(f"Downloaded {filename} to {local_filepath}")
        except Exception as e:
            logging.error(f"Error downloading {filename}: {e}")
    def check_folder_exists_and_not_empty(self, folder_path):
        """Check if a folder exists and is not empty on the remote server."""
        try:
            # Check if folder exists by listing its contents
            file_list = self.sftp_client.listdir(folder_path)
            
            if file_list:
                logging.info(f"The folder '{folder_path}' exists and is not empty.")
                return True  # Folder exists and is not empty
            else:
                logging.info(f"The folder '{folder_path}' exists but is empty.")
                return False  # Folder exists but is empty
        
        except FileNotFoundError:
            logging.error(f"The folder '{folder_path}' does not exist.")
            return False  # Folder does not exist
        
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