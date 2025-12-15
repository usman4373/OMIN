import csv
from pathlib import Path

# ============================================================================
# MAPPINGS AND CONSTANTS
# ============================================================================

# Mapping for force fields
FORCE_FIELDS = {
    '1': {'name': 'CHARMM36', 
          'explicit_files': ['charmm36.xml', 'charmm36/water.xml'], 
          'implicit_compatible': False},
    '2': {'name': 'AMBER14', 
          'explicit_files': ['amber14-all.xml', 'amber14/tip3p.xml'], 
          'main_file': 'amber14-all.xml',
          'implicit_compatible': True},
    '3': {'name': 'AMBER99SB', 
          'explicit_files': ['amber99sb.xml', 'tip3p.xml'], 
          'main_file': 'amber99sb.xml',
          'implicit_compatible': True},
    '4': {'name': 'AMBER03', 
          'explicit_files': ['amber03.xml', 'tip3p.xml'], 
          'main_file': 'amber03.xml',
          'implicit_compatible': True},
    '5': {'name': 'AMBER10', 
          'explicit_files': ['amber10.xml', 'tip3p.xml'], 
          'main_file': 'amber10.xml',
          'implicit_compatible': True}
}

# Mapping for solvent models
SOLVENT_MODELS = {
    '1': {'name': 'Explicit TIP3P', 'type': 'explicit', 'model': 'tip3p'},
    '2': {'name': 'Explicit SPC/E', 'type': 'explicit', 'model': 'spce'},
    '3': {'name': 'Implicit GBn2', 'type': 'implicit', 'model': 'GBn2'},
    '4': {'name': 'Implicit OBC2', 'type': 'implicit', 'model': 'OBC2'}
}

# Force fields compatible with implicit solvent (for auto-switching)
IMPLICIT_COMPATIBLE_FF = ['2', '3', '4', '5']  # AMBER force fields

def check_forcefield_compatibility(selected_ff, selected_sm):
    """Check if selected force field is compatible with solvent model"""
    if selected_sm['type'] == 'implicit' and not selected_ff['implicit_compatible']:
        # Find first compatible force field
        for ff_key in IMPLICIT_COMPATIBLE_FF:
            compatible_ff = FORCE_FIELDS[ff_key]
            print(f"\n⚠️ WARNING: {selected_ff['name']} does not support implicit solvent models.")
            print(f"   Switching to {compatible_ff['name']} force field (compatible with {selected_sm['name']})")
            return compatible_ff
    return selected_ff

def remove_solvent_from_pdb(input_pdb_path, output_pdb_path):
    """Remove solvent (water and ions) from PDB file"""
    solvent_residues = ['HOH', 'WAT', 'TIP', 'SOL', 'NA', 'CL', 'K', 'CA', 'MG']
    
    with open(input_pdb_path, 'r') as f:
        lines = f.readlines()
    
    # Filter out solvent lines
    clean_lines = []
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            resname = line[17:20].strip()
            if resname not in solvent_residues:
                clean_lines.append(line)
        elif line.startswith('TER'):
            clean_lines.append(line)
        elif not line.startswith(('ATOM', 'HETATM', 'TER')):
            # Keep header, title, etc.
            clean_lines.append(line)
    
    # Write clean PDB
    with open(output_pdb_path, 'w') as f:
        f.writelines(clean_lines)

def write_combined_rmsd_csv(rmsd_data, output_dir):
    """
    Write combined RMSD data for all proteins in the requested format
    """
    csv_file = output_dir / "per_residue_rmsd_combined.csv"
    print(f"\nWriting combined RMSD data to: {csv_file}")
    
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Write header with protein names
        header = []
        subheader = []
        
        for protein_name, data in rmsd_data.items():
            header.extend([protein_name, "", ""])
            subheader.extend(["Residue_Name", "Residue_Position", "RMSD_A"])
        
        writer.writerow(header)
        writer.writerow(subheader)
        
        # Find maximum number of residues among all proteins
        max_residues = 0
        for protein_name, data in rmsd_data.items():
            global_rmsd, rmsd_dict, residue_info = data
            max_residues = max(max_residues, len(rmsd_dict))
        
        # Write data row by row
        for i in range(max_residues):
            row = []
            for protein_name, data in rmsd_data.items():
                global_rmsd, rmsd_dict, residue_info = data
                
                # Get sorted residue positions
                sorted_positions = sorted(rmsd_dict.keys())
                
                if i < len(sorted_positions):
                    res_pos = sorted_positions[i]
                    res_name = residue_info.get(res_pos, "UNK")
                    rmsd = rmsd_dict[res_pos]
                    row.extend([res_name, res_pos, f"{rmsd:.4f}"])
                else:
                    # No more residues for this protein
                    row.extend(["", "", ""])
            
            writer.writerow(row)
    
    # Also write a summary CSV with global RMSD values
    summary_file = output_dir / "global_rmsd_summary.csv"
    with open(summary_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Protein", "Global_RMSD_A"])
        for protein_name, data in rmsd_data.items():
            global_rmsd, _, _ = data
            writer.writerow([protein_name, f"{global_rmsd:.4f}"])
    
    print(f"Global RMSD summary saved to: {summary_file}")
    return csv_file, summary_file

def save_metadata(output_dir, params, all_energy_data):
    """Save simulation parameters to a metadata file"""
    metadata_file = output_dir / "simulation_parameters.txt"
    with open(metadata_file, 'w') as f:
        f.write(f"OpenMM Energy Minimization and RMSD Analysis Parameters\n")
        f.write(f"="*60 + "\n")
        f.write(f"Force field: {params['forcefield_info']['name']}\n")
        f.write(f"Solvent model: {params['solvent_info']['name']}\n")
        f.write(f"Minimization iterations: {params['iterations']}\n")
        f.write(f"CPU threads used: {params['threads'] if params['threads'] > 0 else 'All available'}\n")
        
        solvent_type = params['solvent_info']['type']
        if solvent_type == 'explicit':
            f.write(f"Non-bonded method: PME\n")
            f.write(f"Non-bonded cutoff: 1.0 nm\n")
            f.write(f"Solvent padding: 1.0 nm\n")
        else:
            f.write(f"Non-bonded method: NoCutoff\n")
            f.write(f"Implicit solvent model: {params['solvent_info']['model']}\n")
            f.write(f"Solute dielectric: 1.0\n")
            f.write(f"Solvent dielectric: 78.5\n")
        
        f.write(f"pH for adding hydrogens: 7.4\n")
        f.write(f"RMSD analysis performed: {'Yes' if params['run_rmsd'] else 'No'}\n")
        f.write(f"Total proteins processed: {len(all_energy_data)}\n")

        return metadata_file

def validate_directories(input_dir, output_dir):
    """Validate input and output directories"""
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    
    errors = []
    
    # Check input directory
    if not input_path.exists():
        errors.append(f"Input directory does not exist: {input_dir}")
    elif not input_path.is_dir():
        errors.append(f"Input path is not a directory: {input_dir}")
    else:
        # Check for PDB files
        pdb_files = list(input_path.glob("*.pdb")) + list(input_path.glob("*.PDB"))
        if not pdb_files:
            errors.append(f"No PDB files found in input directory: {input_dir}")
    
    # Check output directory (try to create if it doesn't exist)
    try:
        output_path.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        errors.append(f"Cannot create output directory: {output_dir}. Error: {str(e)}")
    
    return errors, input_path, output_path
