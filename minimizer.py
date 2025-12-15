import csv
from pdbfixer import PDBFixer
from openmm.app import PDBFile, ForceField, Simulation, Modeller
from openmm import LangevinIntegrator, Platform
from openmm.unit import kelvin, picosecond, picoseconds, nanometer
import openmm
from openmm.app import PME
import traceback

def process_protein(input_file, output_base_dir, forcefield_info, solvent_info, 
                   iterations, threads, use_gpu, gpu_platform_name, progress_callback=None):
    """Process a single protein file"""
    protein_name = input_file.stem  # Get filename without extension
    
    if progress_callback:
        progress_callback(f"Processing: {protein_name}", 10)
    
    # Step 1: Fix PDB using PDBFixer (without saving intermediate file)
    if progress_callback:
        progress_callback("Fixing PDB...", 20)
    
    fixer = PDBFixer(filename=str(input_file))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    
    # Step 2: Prepare system (directly from fixer without saving)
    if progress_callback:
        progress_callback("Setting up system...", 30)
    
    topology = fixer.topology
    positions = fixer.positions
    
    modeller = Modeller(topology, positions)
    
    # Load force field based on solvent type
    try:
        if solvent_info['type'] == 'explicit':
            # For explicit solvent, use the explicit_files
            forcefield_files = forcefield_info['explicit_files']
            
            # If SPC/E is selected, we need to use the correct water model file
            if solvent_info['model'] == 'spce':
                # Check if we need to adjust the water model file
                if forcefield_info['name'] == 'CHARMM36':
                    forcefield_files = ['charmm36.xml', 'charmm36/spce.xml']
                elif forcefield_info['name'] == 'AMBER14':
                    forcefield_files = ['amber14-all.xml', 'amber14/spce.xml']
                else:
                    # For other AMBER force fields, spce.xml is in the root
                    forcefield_files = [forcefield_info['explicit_files'][0], 'spce.xml']
            
            forcefield = ForceField(*forcefield_files)
        else:
            # For implicit solvent, use main force field file + appropriate implicit solvent file
            main_file = forcefield_info.get('main_file', forcefield_info['explicit_files'][0])
            
            # Determine which implicit solvent file to use
            if solvent_info['model'] == 'GBn2':
                implicit_file = 'implicit/gbn2.xml'
            elif solvent_info['model'] == 'OBC2':
                implicit_file = 'implicit/obc2.xml'
            else:
                # Default to GBn2
                implicit_file = 'implicit/gbn2.xml'
            
            forcefield_files = [main_file, implicit_file]
            print(f"  Loading force field files: {forcefield_files}")
            forcefield = ForceField(*forcefield_files)
    except Exception as e:
        print(f"Error loading force field: {e}")
        # Fallback: try with just the main force field file
        try:
            if solvent_info['type'] == 'explicit':
                forcefield = ForceField(forcefield_info['explicit_files'][0])
            else:
                forcefield = ForceField(forcefield_info['explicit_files'][0])
                print("Warning: Using force field without specific implicit solvent parameters")
        except Exception as e2:
            print(f"Fallback also failed: {e2}")
            raise
    
    # Handle different solvent types
    if solvent_info['type'] == 'explicit':
        if progress_callback:
            progress_callback(f"Adding explicit solvent ({solvent_info['model']})...", 40)
        
        modeller.addSolvent(forcefield, model=solvent_info['model'], padding=1.0*nanometer)
        
        # Create system with explicit solvent (PME for electrostatics)
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=PME,
            nonbondedCutoff=1.0*nanometer
        )
    else:
        # Implicit solvent
        if progress_callback:
            progress_callback(f"Setting up implicit solvent ({solvent_info['model']})...", 40)
        
        # Create system with implicit solvent
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=openmm.app.NoCutoff,
            soluteDielectric=1.0,
            solventDielectric=78.5
        )
    
    # Step 3: Create simulation
    if progress_callback:
        progress_callback("Setting up simulation...", 50)
    
    # Create integrator
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

    # Create simulation with platform properties
    if use_gpu and gpu_platform_name:
        try:
            platform = Platform.getPlatformByName(gpu_platform_name)
            # Use default GPU properties (OpenMM will choose the fastest available GPU)
            properties = {}  
            simulation = Simulation(modeller.topology, system, integrator, platform, properties)
            print(f"  Using {gpu_platform_name} GPU acceleration")
        except Exception as e:
            print(f"  GPU setup failed: {e}. Falling back to CPU.")
            if threads > 0:
                platform = Platform.getPlatformByName('CPU')
                properties = {'Threads': str(threads)}
                simulation = Simulation(modeller.topology, system, integrator, platform, properties)
            else:
                simulation = Simulation(modeller.topology, system, integrator)
    elif threads > 0:
        platform = Platform.getPlatformByName('CPU')
        properties = {'Threads': str(threads)}
        simulation = Simulation(modeller.topology, system, integrator, platform, properties)
        print(f"  Using CPU with {threads} thread(s)")
    else:
        simulation = Simulation(modeller.topology, system, integrator)
        print("  Using CPU with default threads")
    
    simulation.context.setPositions(modeller.positions)
    
    # Step 4: Get initial energy
    if progress_callback:
        progress_callback("Calculating initial energy...", 60)
    
    state = simulation.context.getState(getEnergy=True)
    initial_energy = state.getPotentialEnergy().value_in_unit(openmm.unit.kilojoule_per_mole)
    
    # Step 5: Minimize
    if progress_callback:
        progress_callback(f"Minimizing energy ({iterations} iterations)...", 70)
    
    simulation.minimizeEnergy(maxIterations=iterations)
    
    # Step 6: Get minimized energy
    if progress_callback:
        progress_callback("Getting minimized energy...", 80)
    
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    minimized_energy = state.getPotentialEnergy().value_in_unit(openmm.unit.kilojoule_per_mole)
    minimized_positions = state.getPositions()
    
    # Import remove_solvent_from_pdb from utils
    from utils import remove_solvent_from_pdb
    
    # Create output directories
    if solvent_info['type'] == 'explicit':
        solvent_dir = output_base_dir / "minimized-solvent"
        clean_dir = output_base_dir / "minimized-clean"
        
        solvent_dir.mkdir(exist_ok=True)
        clean_dir.mkdir(exist_ok=True)
        
        # Step 7: Save results
        # Save with solvent
        solvent_file = solvent_dir / f"{protein_name}_minimized.pdb"
        print(f"Saving with solvent to: {solvent_file}")
        with open(solvent_file, 'w') as f:
            PDBFile.writeFile(modeller.topology, minimized_positions, f)
        
        # Save without solvent
        clean_file = clean_dir / f"{protein_name}_minimized.pdb"
        print(f"Saving without solvent to: {clean_file}")
        
        # First save the minimized structure to a temporary file
        temp_file = clean_dir / f"{protein_name}_temp.pdb"
        with open(temp_file, 'w') as f:
            PDBFile.writeFile(modeller.topology, minimized_positions, f)
        
        # Remove solvent and save clean version
        remove_solvent_from_pdb(temp_file, clean_file)
        temp_file.unlink()  # Remove temporary file
        
    else:
        # For implicit solvent, just save the minimized structure
        implicit_dir = output_base_dir / "minimized-implicit"
        implicit_dir.mkdir(exist_ok=True)
        
        clean_file = implicit_dir / f"{protein_name}_minimized.pdb"
        print(f"Saving minimized structure to: {clean_file}")
        with open(clean_file, 'w') as f:
            PDBFile.writeFile(modeller.topology, minimized_positions, f)
    
    if progress_callback:
        progress_callback("Minimization complete!", 100)
    
    print(f"Initial energy: {initial_energy:.4f} kJ/mol")
    print(f"Minimized energy: {minimized_energy:.4f} kJ/mol")
    print(f"Energy difference: {initial_energy - minimized_energy:.4f} kJ/mol")
    
    return {
        'filename': protein_name,
        'raw_file': input_file,
        'minimized_file': clean_file,
        'initial_energy': initial_energy,
        'final_energy': minimized_energy,
        'energy_difference': initial_energy - minimized_energy,
        'clean_file': clean_file
    }

def save_energy_data(all_energy_data, output_dir):
    """Save energy data to CSV file"""
    energies_file = output_dir / "energies.csv"
    print(f"\nSaving energy data to: {energies_file}")
    
    with open(energies_file, 'w', newline='') as csvfile:
        # Write header
        writer = csv.writer(csvfile)
        writer.writerow(['filename', 'initial_energy (kJ/mol)', 'final_energy (kJ/mol)', 
                       'energy_difference (kJ/mol)'])
        
        # Write data for each protein
        for data in all_energy_data:
            # Format the values to 4 decimal places if they're numbers
            if data['initial_energy'] != 'ERROR':
                initial_str = f"{data['initial_energy']:.4f}"
                final_str = f"{data['final_energy']:.4f}"
                diff_str = f"{data['energy_difference']:.4f}"
            else:
                initial_str = data['initial_energy']
                final_str = data['final_energy']
                diff_str = data['energy_difference']
            
            writer.writerow([
                data['filename'],
                initial_str,
                final_str,
                diff_str
            ])
    
    return energies_file

def process_all_proteins(pdb_files, output_dir, forcefield_info, solvent_info, 
                        iterations, threads, use_gpu, gpu_platform_name, progress_callback=None):
    """Process all PDB files"""
    all_energy_data = []
    protein_pairs = {}
    
    for i, pdb_file in enumerate(pdb_files, 1):
        if progress_callback:
            progress_callback(f"Processing file {i}/{len(pdb_files)}: {pdb_file.name}", 
                            int((i-1) * 100 / len(pdb_files)))
        
        try:
            result = process_protein(
                pdb_file, 
                output_dir, 
                forcefield_info, 
                solvent_info, 
                iterations, 
                threads,
                use_gpu,
                gpu_platform_name,
                progress_callback
            )
            all_energy_data.append(result)
            
            # Store protein pair for RMSD analysis
            protein_name = pdb_file.stem
            protein_pairs[protein_name] = {
                'raw': pdb_file,
                'minimized': result['clean_file']
            }
            
        except Exception as e:
            print(f"Error processing {pdb_file.name}: {e}")
            traceback.print_exc()
            # Add placeholder data for failed protein
            all_energy_data.append({
                'filename': pdb_file.stem,
                'raw_file': pdb_file,
                'minimized_file': None,
                'initial_energy': 'ERROR',
                'final_energy': 'ERROR',
                'energy_difference': 'ERROR',
                'clean_file': None
            })
    
    return all_energy_data, protein_pairs