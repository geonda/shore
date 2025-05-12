import os
import sys
import numpy as np
from ase.build import bulk
from ase.io import write as ase_write
from pathlib import Path

# Add the parent directory to the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from shore import AtomicStructure
from shore.input_manager import Input, Matter, Light
from shore.remote import RemoteServerManager
from shore.calculation import Calculation

def create_test_structure():
    """Create a slightly distorted Si bulk structure for testing relaxation."""
    # Create a slightly distorted silicon structure
    si = bulk('Si', 'diamond', a=5.43, cubic=True)
    
    # Add some distortion to test relaxation
    cell = si.get_cell()
    cell[0][0] *= 1.02  # 2% distortion in x
    si.set_cell(cell, scale_atoms=True)
    
    # Slightly displace one atom
    positions = si.get_positions()
    positions[0][0] += 0.1  # Displace first atom by 0.1 Å
    si.set_positions(positions)
    
    return si

def setup_test_calculation():
    """Setup a test calculation with the distorted structure."""
    # Create test directories
    test_dir = Path('test_relax')
    test_dir.mkdir(exist_ok=True)
    
    # Save the structure
    atoms = create_test_structure()
    structure_file = test_dir / 'si_distorted.xyz'
    ase_write(structure_file, atoms)
    
    # Create structure object using AtomicStructure
    structure = AtomicStructure(filename=str(structure_file))
    
    # Create Matter object with QE parameters
    matter = Matter(structure=structure, 
                   pseudo_dir='./tests/pseudo',
                   ecutwfc=40,
                   ecutrho=320)
    
    # Create Light object (empty for relaxation)
    light = Light()
    
    # Create Input object
    input_obj = Input(name='si_test_calc', 
                     target='qe',  # Set target to QE for relaxation
                     matter=matter,
                     light=light)
    
    # Create server manager (update with your server details or set to None for local)
    server = None  # RemoteServerManager(hostname='your_host', username='your_username')
    
    # Create calculation object
    calc = Calculation(server=server, input=input_obj)
    
    return calc

def test_relaxation():
    """Run the test case."""
    print("Starting relaxation test...")
    
    # Setup calculation
    calc = setup_test_calculation()
    
    # Setup relaxation
    print("Setting up QE relaxation...")
    calc.setup_relax()
    
    # Print the QE input for verification
    relax_input_file = os.path.join(calc.local_dir, 'relax', 'relax.in')
    print("\nGenerated QE input file contents:")
    print("-" * 50)
    with open(relax_input_file, 'r') as f:
        print(f.read())
    print("-" * 50)
    
    # Run relaxation if QE is available
    try:
        print("\nAttempting to run relaxation...")
        calc.run_relax(monitor=True)
        
        # Check results if relaxation completed
        relaxed_structure_file = os.path.join(calc.local_dir, 'relax', 'relaxed.xyz')
        if os.path.exists(relaxed_structure_file):
            print("\nRelaxation completed successfully!")
            print("Relaxed structure saved to:", relaxed_structure_file)
        else:
            print("\nRelaxation output not found. Check if QE ran successfully.")
            
    except Exception as e:
        print("\nError running relaxation:", str(e))
        print("Make sure Quantum ESPRESSO is installed and accessible.")

if __name__ == '__main__':
    test_relaxation() 