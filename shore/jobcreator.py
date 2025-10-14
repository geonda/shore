class JobScriptCreator:
    def __init__(self, job_name='test.job', time_limit='01:00:00', ncores=1):
        self.job_name = job_name
        self.output_file = 'out'
        self.error_file = 'err'
        self.time_limit = time_limit
        self.ntasks_per_node=ncores
        self.cpus_per_task=1
        self.nodes=1

    def generate_script(self,path='./', command=None, module='', extra=''):
        """Generate the job.sh script with the specified executable."""
        script_content = [
            "#!/bin/bash",
            f"#SBATCH --job-name={self.job_name}",
            f"#SBATCH --output={self.output_file}",
            f"#SBATCH --error={self.error_file}",
            f"#SBATCH --ntasks-per-node={self.ntasks_per_node}",
            f"#SBATCH --cpus-per-task={self.cpus_per_task}",
            f"#SBATCH --nodes={self.nodes}",
            f"#SBATCH --time={self.time_limit}",
            f"{extra}\n",  # Command to execute the program
            f"{module}\n",  # Command to execute the program
            f"{command}\n"  # Command to execute the program
        ]

        with open(f"{path}/job.sh", "w") as f:
            f.write("\n".join(script_content))
        
        # print(f"Job script 'job.sh' created successfully.")
