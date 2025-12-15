import pymol
from pymol import cmd

def get_all_residues_with_info(object_name):
    """Get all residues with their information"""
    residues = []
    cmd.iterate(f"{object_name} and name CA", 
                "residues.append((chain, resi, resn))", 
                space={'residues': residues})
    return residues

def compute_per_residue_rmsd_fixed(raw_obj, min_obj):
    """
    Compute per-residue RMSD - fixed version that handles chain differences
    """
    rmsd_dict = {}
    residue_info = {}  # Store residue info for CSV
    
    # Get all residues from both structures
    raw_residues = get_all_residues_with_info(raw_obj)
    min_residues = get_all_residues_with_info(min_obj)
    
    # Create a mapping: residue number -> (chain, resn) for each structure
    raw_map = {resi: (chain, resn) for chain, resi, resn in raw_residues}
    min_map = {resi: (chain, resn) for chain, resi, resn in min_residues}
    
    # Find common residue numbers
    common_residues = set(raw_map.keys()) & set(min_map.keys())
    common_residues = sorted(common_residues, key=lambda x: int(x))
    
    print(f"\nFound {len(common_residues)} common residues")
    
    for resi in common_residues:
        try:
            raw_chain, raw_resn = raw_map[resi]
            min_chain, min_resn = min_map[resi]
            
            # Store residue info for CSV
            residue_info[int(resi)] = raw_resn
            
            # Check if residue names match
            if raw_resn != min_resn:
                print(f"Warning: Residue {resi} has different names: {raw_resn} vs {min_resn}")
            
            # Create selections with specific chains
            raw_sel = f"{raw_obj} and chain {raw_chain} and resi {resi}"
            min_sel = f"{min_obj} and chain {min_chain} and resi {resi}"
            
            # Check if selections have atoms
            n_raw = cmd.count_atoms(raw_sel)
            n_min = cmd.count_atoms(min_sel)
            
            if n_raw == 0 or n_min == 0:
                print(f"Warning: No atoms found for residue {resi}")
                continue
            
            # Use align with cycles=0 to get RMSD without further fitting
            # This matches atoms by name
            rmsd_result = cmd.align(min_sel, raw_sel, cycles=0)
            
            if rmsd_result and len(rmsd_result) > 0:
                rmsd = rmsd_result[0]
                rmsd_dict[int(resi)] = rmsd
                
                if int(resi) % 10 == 0:
                    print(f"  Residue {resi} ({raw_resn}): {rmsd:.3f} Å")
            else:
                print(f"Warning: Could not compute RMSD for residue {resi}")
                
        except Exception as e:
            print(f"Error for residue {resi}: {str(e)}")
            continue
    
    return rmsd_dict, residue_info

def compute_rmsd_analysis(raw_pdb_path, min_pdb_path, protein_name):
    """
    Compute RMSD analysis for a single protein pair
    Returns: global_rmsd, per_residue_rmsd_dict, residue_info_dict
    """
    # Initialize PyMOL if not already done
    try:
        # Check if PyMOL is already launched
        cmd.get_names()
    except:
        # Launch PyMOL WITHOUT exiting flags
        pymol.finish_launching(['pymol', '-cqi'])  # -i keeps it interactive and prevents auto-exit
    
    # Load structures with unique names
    raw_obj = f"raw_{protein_name}"
    min_obj = f"min_{protein_name}"
    
    cmd.load(str(raw_pdb_path), raw_obj)
    cmd.load(str(min_pdb_path), min_obj)
    
    # Align globally
    print(f"  Aligning structures for {protein_name}...")
    global_rmsd = cmd.align(min_obj, raw_obj)[0]
    print(f"  Global RMSD: {global_rmsd:.3f} Å")
    
    # Compute per-residue RMSD
    print(f"  Computing per-residue RMSD for {protein_name}...")
    rmsd_dict, residue_info = compute_per_residue_rmsd_fixed(raw_obj, min_obj)
    
    if not rmsd_dict:
        print(f"  WARNING: No per-residue RMSD values computed for {protein_name}!")
        print(f"  Using global RMSD for all residues...")
        # Create dummy data with global RMSD for all residues
        raw_residues = get_all_residues_with_info(raw_obj)
        rmsd_dict = {}
        for chain, resi, resn in raw_residues:
            rmsd_dict[int(resi)] = global_rmsd
            residue_info[int(resi)] = resn
    
    # Clean up objects
    cmd.delete(raw_obj)
    cmd.delete(min_obj)
    
    return global_rmsd, rmsd_dict, residue_info

def create_custom_color_gradient():
    """
    Create a custom color gradient from greencyan to red
    """
    # Define color stops for our gradient
    
    # First, ensure we have the base colors defined
    # Greencyan is like a teal color (cyan + green)
    cmd.set_color("greencyan", [0.0, 0.7, 0.7])  # RGB: teal/cyan-green
    
    # Create intermediate colors for smooth gradient
    # From greencyan (0, 0.7, 0.7) to orange (1, 0.5, 0) to red (1, 0, 0)
    
    # Color at 0.25: between greencyan and orange
    cmd.set_color("gradient_01", [0.5, 0.6, 0.35])
    # Color at 0.5: orange
    cmd.set_color("gradient_02", [1.0, 0.5, 0.0])
    # Color at 0.75: between orange and red
    cmd.set_color("gradient_03", [1.0, 0.25, 0.0])
    
    return ["greencyan", "gradient_01", "gradient_02", "gradient_03", "red"]

def color_by_rmsd_custom(rmsd_dict, min_copy_obj):
    """
    Apply custom coloring based on RMSD values
    """
    if not rmsd_dict or len(rmsd_dict) <= 1:
        return
    
    values = list(rmsd_dict.values())
    min_val = min(values)
    max_val = max(values)
    
    print(f"  RMSD range: {min_val:.3f} to {max_val:.3f} Å")
    

    # Create custom color gradient
    color_names = create_custom_color_gradient()
    
    # Normalize RMSD values to 0-1 range for coloring
    if max_val > 0:
        if max_val - min_val > 0.001:  # Avoid division by zero
            # Create color stops at appropriate positions
            # Map RMSD values to color positions
            stops = [0.00, 0.01, 0.10, 0.20, 0.30]
            
            # Use ramp_new to create custom gradient
            # But we need to avoid displaying it, so we'll hide it immediately
            cmd.ramp_new("custom_gradient", min_copy_obj, stops, color_names)
            
            # Hide the ramp (colorbar) from view
            cmd.hide("everything", "custom_gradient")
            
            # Alternative: Use spectrum with custom palette
            # Create a string of color names for spectrum
            color_string = " ".join(color_names)
            cmd.spectrum("b", color_string, min_copy_obj, 
                        minimum=0.00, maximum=0.30)
        else:
            # All values are the same, use greencyan
            cmd.color("greencyan", min_copy_obj)
    else:
        # All RMSD values are 0, use greencyan
        cmd.color("greencyan", min_copy_obj)

def visualize_rmsd_analysis(raw_pdb_path, min_pdb_path, protein_name, output_dir, image_width=1080, image_height=1080, image_dpi=300):
    """
    Create visualization for RMSD analysis
    """
    print(f"\nCreating visualization for {protein_name}...")
    
    # Load structures
    raw_obj = f"raw_{protein_name}"
    min_obj = f"min_{protein_name}"
    
    cmd.load(str(raw_pdb_path), raw_obj)
    cmd.load(str(min_pdb_path), min_obj)
    
    # Align globally
    global_rmsd = cmd.align(min_obj, raw_obj)[0]
    
    # Create copy for coloring
    min_copy_obj = f"min_copy_{protein_name}"
    cmd.create(min_copy_obj, min_obj)
    
    # Compute per-residue RMSD
    rmsd_dict, residue_info = compute_per_residue_rmsd_fixed(raw_obj, min_obj)
    
    if not rmsd_dict:
        print(f"  WARNING: No per-residue RMSD values computed!")
        # Set all B-factors to global RMSD
        cmd.alter(min_copy_obj, f"b = {global_rmsd}")
        # Create dummy data for summary
        rmsd_dict = {1: global_rmsd}
        # Get residue info from raw structure
        raw_residues = get_all_residues_with_info(raw_obj)
        for chain, resi, resn in raw_residues:
            residue_info[int(resi)] = resn
    else:
        print(f"  Computed RMSD for {len(rmsd_dict)} residues")
        # Set B-factors for all atoms based on RMSD
        for resi, rmsd in rmsd_dict.items():
            # Find chain for this residue
            raw_residues = get_all_residues_with_info(raw_obj)
            for chain, res, resn in raw_residues:
                if int(res) == resi:
                    cmd.alter(f"{min_copy_obj} and chain {chain} and resi {resi}", f"b = {rmsd}")
                    break
    
    # Step 1: Color raw structure with greencyan
    cmd.set_color("greencyan", [0.0, 0.7, 0.7])  # RGB: teal/cyan-green
    cmd.color("greencyan", raw_obj)
    
    # Step 2: Color minimized structure with custom gradient based on RMSD
    color_by_rmsd_custom(rmsd_dict, min_copy_obj)
    
    # Visualization settings
    cmd.hide("everything")
    cmd.show("cartoon", raw_obj)
    cmd.show("cartoon", min_copy_obj)
    cmd.set("cartoon_transparency", 0.5, raw_obj)
    cmd.set("cartoon_transparency", 0.0, min_copy_obj)
    cmd.set('cartoon_highlight_color', 'grey35')
    cmd.set("cartoon_fancy_helices", 1)
    
    # Rendering parameters
    cmd.set('ray_trace_mode', 3)
    cmd.set('ray_trace_gain', 0.6)
    cmd.set('ambient', 0.3)
    cmd.set('direct', 0.6)
    cmd.set('specular', 0.2)
    cmd.set('shininess', 200)
    cmd.set('ray_shadow', 1)
    cmd.set("light_count", 2)
    cmd.set('antialias', 4)
    cmd.set("depth_cue", 0)
    cmd.set("fog", 0)
    
    # Orient and set background
    cmd.orient()
    cmd.zoom()
    cmd.bg_color("white")
    
    # Remove any existing ramps/colorbars
    cmd.delete("rmsd_ramp")
    cmd.delete("custom_gradient")
    
    # Create output directory for visualizations
    viz_dir = output_dir / "visualizations"
    viz_dir.mkdir(exist_ok=True)
    
    # Save session
    session_file = viz_dir / f"{protein_name}.pse"
    cmd.save(str(session_file))
    
    # Save image
    image_file = viz_dir / f"{protein_name}.png"
    cmd.ray(image_width, image_height)
    cmd.png(str(image_file), width=image_width, height=image_height, dpi=image_dpi, ray=1)
    
    print(f"  Visualization saved to: {viz_dir}/")
    print(f"    - {protein_name}.pse (PyMOL session)")
    print(f"    - {protein_name}.png (high-resolution image)")
    print(f"  Visualization complete for {protein_name}")
    
    # Clean up objects but DON'T exit PyMOL
    cmd.delete(raw_obj)
    cmd.delete(min_obj)
    cmd.delete(min_copy_obj)
    
    return global_rmsd, rmsd_dict, residue_info

def run_rmsd_analysis(protein_pairs, output_dir, progress_callback=None, image_width=1080, image_height=1080, image_dpi=300):
    """
    Run RMSD analysis for all protein pairs
    """
    # Initialize PyMOL once - BUT don't let it exit
    print("Initializing PyMOL for visualization...")
    try:
        # Check if PyMOL is already running
        cmd.get_names()
        print("PyMOL already initialized")
    except:
        # Launch PyMOL with -i flag to prevent auto-exit
        pymol.finish_launching(['pymol', '-cqi'])
        print("PyMOL launched in interactive mode (won't auto-exit)")
    
    # Dictionary to store RMSD data for all proteins
    all_rmsd_data = {}
    
    # Process each protein pair for RMSD analysis
    successful_rmsd = 0
    total_pairs = len(protein_pairs)
    
    for idx, (protein_name, paths) in enumerate(protein_pairs.items(), 1):
        if progress_callback:
            progress_callback(f"Analyzing RMSD for: {protein_name}", 
                            int((idx-1) * 100 / total_pairs))
        
        print(f"\nAnalyzing RMSD for: {protein_name}")
        print("-"*40)
        
        if paths['minimized'] and paths['minimized'].exists():
            try:
                # Compute RMSD analysis
                global_rmsd, rmsd_dict, residue_info = compute_rmsd_analysis(
                    paths['raw'], 
                    paths['minimized'], 
                    protein_name
                )
                
                # Store data
                all_rmsd_data[protein_name] = (global_rmsd, rmsd_dict, residue_info)
                
                # Create visualization
                viz_global_rmsd, viz_rmsd_dict, viz_residue_info = visualize_rmsd_analysis(
                    paths['raw'],
                    paths['minimized'],
                    protein_name,
                    output_dir,
                    image_width,
                    image_height,
                    image_dpi
                )
                
                successful_rmsd += 1
                
                # Print summary for this protein
                if rmsd_dict:
                    values = list(rmsd_dict.values())
                    avg_rmsd = sum(values) / len(values)
                    max_res = max(rmsd_dict.items(), key=lambda x: x[1])
                    min_res = min(rmsd_dict.items(), key=lambda x: x[1])
                    
                    print(f"\n  RMSD Summary for {protein_name}:")
                    print(f"    Global RMSD:          {global_rmsd:.3f} Å")
                    print(f"    Average per-residue:  {avg_rmsd:.3f} Å")
                    print(f"    Maximum RMSD:         {max_res[1]:.3f} Å (Residue {max_res[0]})")
                    print(f"    Minimum RMSD:         {min_res[1]:.3f} Å (Residue {min_res[0]})")
                
            except Exception as e:
                print(f"  Error in RMSD analysis for {protein_name}: {e}")
                import traceback
                traceback.print_exc()
        else:
            print(f"  Skipping {protein_name}: Minimized file not found")
    
    if progress_callback:
        progress_callback("RMSD analysis complete!", 100)
    
    return all_rmsd_data, successful_rmsd