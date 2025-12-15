import streamlit as st
import os
from pathlib import Path
import traceback
import base64

# Import modules
from utils import FORCE_FIELDS, SOLVENT_MODELS, check_forcefield_compatibility, save_metadata
from minimizer import process_all_proteins, save_energy_data
from visualization import run_rmsd_analysis

# Set page configuration
st.set_page_config(
    page_title="OMIN",
    page_icon="⚡︎",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Initialize session state
if 'processing' not in st.session_state:
    st.session_state.processing = False
if 'results' not in st.session_state:
    st.session_state.results = None
if 'progress' not in st.session_state:
    st.session_state.progress = 0
if 'progress_text' not in st.session_state:
    st.session_state.progress_text = ""
if 'processing_stopped' not in st.session_state:
    st.session_state.processing_stopped = False
if 'processed_files' not in st.session_state:
    st.session_state.processed_files = []
if 'visualization_images' not in st.session_state:
    st.session_state.visualization_images = []
if 'current_progress' not in st.session_state:
    st.session_state.current_progress = []
if 'file_status' not in st.session_state:
    st.session_state.file_status = {}
if 'current_file' not in st.session_state:
    st.session_state.current_file = None

def main():
    """Display main header"""

# Add icon
with open("icon/icon.png", "rb") as f:
    data = f.read()
    data_base64 = base64.b64encode(data).decode()

    # app name
    st.markdown(
    f"""
    <div style='display: flex; align-items: center;'>
        <h1 style='margin-bottom: 0; font-size: 80px; font-weight: 800;
        '>OMIN</h1>
        <img src='data:image/png;base64,{data_base64}' width='150' style='margin-bottom: 16px;'>
    </div>
    """,
    unsafe_allow_html=True
)

    st.markdown(
    """
    <h5 style="font-weight: 300;">
        <span style="font-weight:700; color: var(--primary-color);">O</span>penMM
        <span style="font-weight:700; color: var(--secondary-color);">MIN</span>imizer
        for Protein Energy Minimizations
    </h5>
    """,
    unsafe_allow_html=True
)

    st.markdown("## ⚙️ Simulation Parameters")

    # Create columns for better organization
    col1, col2 = st.columns(2)

    with col1:
        # Input directory path
        input_dir = st.text_input(
            "**Input Directory Path**",
            placeholder="Add/input/directory/path",
            help="Directory containing .pdb files to process"
        )

        # Force field selection
        forcefield_options = {ff['name']: key for key, ff in FORCE_FIELDS.items()}
        selected_ff_name = st.selectbox(
            "**Select Force Field:**",
            options=list(forcefield_options.keys()),
            index=1
        )
        selected_ff_key = forcefield_options[selected_ff_name]
        selected_ff = FORCE_FIELDS[selected_ff_key]

        # Minimization parameters
        iterations = st.number_input(
            "**Number of Iterations:**",
            min_value=10,
            max_value=10000,
            value=500,
            step=100
        )

    with col2:
        # Output directory path
        output_dir = st.text_input(
            "**Output Directory Path:**",
            placeholder="Add/output/directory/path",
            help="Directory where results will be saved"
        )

        # Solvent model selection
        solvent_options = {sm['name']: key for key, sm in SOLVENT_MODELS.items()}
        selected_sm_name = st.selectbox(
            "**Select Solvent Model:**",
            options=list(solvent_options.keys()),
            index=2
        )
        selected_sm_key = solvent_options[selected_sm_name]
        selected_sm = SOLVENT_MODELS[selected_sm_key]

        # Hardware selection
        hardware_choice = st.radio(
            "**Hardware Acceleration:**",
            ["CPU", "GPU"],
            index=0,
            horizontal=True
        )
        use_gpu = (hardware_choice == "GPU")

        if not use_gpu:
            threads = st.number_input(
                "**CPU threads** (0 for all available):",
                min_value=0,
                max_value=os.cpu_count() or 128,
                value=0,
                step=2
            )
        else:
            threads = 0

    # Check compatibility
    selected_ff = check_forcefield_compatibility(selected_ff, selected_sm)
    if selected_ff['name'] != selected_ff_name:
        st.info(f"⚠️ Auto-switched to {selected_ff['name']} for compatibility with {selected_sm['name']}")

    # Visualizations section
    st.markdown("### 🎨 Visualizations")
    col_viz1, col_viz2 = st.columns([1, 2])

    with col_viz1:
        run_visualizations = st.checkbox("Run visualizations after minimization", value=True)

    with col_viz2:
        if run_visualizations:
            viz_col1, viz_col2, viz_col3 = st.columns(3)
            with viz_col1:
                image_width = st.number_input("Image Width", min_value=100, max_value=4000, value=1080, step=500)
            with viz_col2:
                image_height = st.number_input("Image Height", min_value=100, max_value=4000, value=1080, step=500)
            with viz_col3:
                image_dpi = st.number_input("Image DPI", min_value=72, max_value=1200, value=300, step=100)

    # Run button at the bottom of parameters section
    st.markdown("---")
    # Check if input directory exists and has PDB files
    input_path = Path(input_dir)
    pdb_files_exist = False
    if input_path.exists() and input_path.is_dir():
        pdb_files = list(input_path.glob("*.pdb")) + list(input_path.glob("*.PDB"))
        pdb_files_exist = len(pdb_files) > 0

    run_button = st.button(
        "🚀 Run Energy Minimizer",
        type="primary",
        disabled=st.session_state.processing or not pdb_files_exist
    )
    
    # Run processing when button is clicked
    if run_button and not st.session_state.processing:
        # Validate input and output directories
        input_path = Path(input_dir)
        output_path = Path(output_dir)
        
        if not input_path.exists() or not input_path.is_dir():
            st.error(f"❌ Input directory does not exist: {input_dir}")
            st.stop()
        
        pdb_files = list(input_path.glob("*.pdb")) + list(input_path.glob("*.PDB"))
        if not pdb_files:
            st.error(f"❌ No PDB files found in {input_dir}")
            st.stop()
        
        # Create output directory if it doesn't exist
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Set processing state
        st.session_state.processing = True
        st.session_state.progress = 0
        st.session_state.progress_text = "Starting..."
        
        st.markdown("## ⏳ Processing")
        
        # Initialize file status tracking
        pdb_filenames = [f.stem for f in pdb_files]
        if 'file_status' not in st.session_state:
            st.session_state.file_status = {}
            # Initialize all files with empty status dictionaries
            for name in pdb_filenames:
                st.session_state.file_status[name] = {'minimization': 'In Queue', 'visualization': ''}
        if 'current_file' not in st.session_state:
            st.session_state.current_file = None
        
        # Create a placeholder for status display
        status_container = st.empty()
        
        def update_progress(text, value):
            """Update progress and file status for both minimization and visualization"""
            # Store progress values
            st.session_state.progress = value
            st.session_state.progress_text = text
            
            # Debug: print the text to see what messages we're getting
            print(f"Progress update: '{text}'")
            
            # Parse the text to update file status
            # Phase 1: Minimization
            if text.startswith("Processing file"):
                # Extract filename from message like "Processing file X/Y: filename.pdb"
                try:
                    filename = text.split(": ")[1].replace(".pdb", "").replace(".PDB", "")
                    if filename not in st.session_state.file_status:
                        st.session_state.file_status[filename] = {'minimization': 'In Queue', 'visualization': ''}
                    st.session_state.current_file = filename
                    # Mark current file as processing for minimization
                    st.session_state.file_status[filename]['minimization'] = "processing"
                    # Mark previous file as complete if it exists
                    for f in pdb_filenames:
                        if (f != filename and 
                            f in st.session_state.file_status and
                            st.session_state.file_status[f].get('minimization') == "processing"):
                            st.session_state.file_status[f]['minimization'] = "✅"
                    print(f"  -> Marked {filename} as processing (minimization)")
                except Exception as e:
                    print(f"Error parsing filename: {e}")
            
            elif "Minimization complete!" in text or "Minimization complete" in text:
                # Mark current file as complete for minimization
                if st.session_state.current_file and st.session_state.current_file in st.session_state.file_status:
                    st.session_state.file_status[st.session_state.current_file]['minimization'] = "✅"
                    print(f"  -> Marked {st.session_state.current_file} as ✅ (minimization complete)")
            
            # Phase 2: Visualization - Look for start of visualization
            elif "Creating visualization for" in text or "Analyzing RMSD for" in text:
                # Extract protein name from message
                try:
                    if "Creating visualization for" in text:
                        protein_name = text.replace("Creating visualization for ", "").replace("...", "").strip()
                    elif "Analyzing RMSD for" in text:
                        protein_name = text.replace("Analyzing RMSD for: ", "").strip()
                    
                    if protein_name and protein_name not in st.session_state.file_status:
                        st.session_state.file_status[protein_name] = {'minimization': '✅', 'visualization': 'In Queue'}
                    
                    if protein_name:
                        st.session_state.current_file = protein_name
                        # Mark as processing for visualization
                        st.session_state.file_status[protein_name]['visualization'] = "processing"
                        # Mark previous file as complete if it exists (same logic as minimization)
                        for f in pdb_filenames:
                            if (f != protein_name and 
                                f in st.session_state.file_status and
                                st.session_state.file_status[f].get('visualization') == "processing"):
                                st.session_state.file_status[f]['visualization'] = "✅"
                        print(f"  -> Marked {protein_name} as processing (visualization)")
                except Exception as e:
                    print(f"Error parsing visualization message: {e}")
            
            # Look for visualization completion messages - handle the last file
            elif (
                "Visualization saved to:" in text
                or "RMSD Summary for" in text
                or "Visualization complete for" in text
            ):
                # Mark current file as complete for visualization
                if (
                    st.session_state.current_file
                    and st.session_state.current_file in st.session_state.file_status
                ):
                    st.session_state.file_status[
                        st.session_state.current_file
                    ]['visualization'] = "✅"
                    print(
                        f"  -> Marked {st.session_state.current_file} as ✅ (visualization complete)"
                    )
            
            elif "Processing complete!" in text:
                if (
                    st.session_state.current_file
                    and st.session_state.current_file in st.session_state.file_status
                    and st.session_state.file_status[
                        st.session_state.current_file
                    ].get("visualization") == "processing"
                ):
                    st.session_state.file_status[
                        st.session_state.current_file
                    ]['visualization'] = "✅"
                    print(
                        f"  -> Finalized {st.session_state.current_file} as ✅ (last visualization)"
                    )

            
            # When visualization phase starts, mark all successfully minimized files as "In Queue" for visualization
            elif "Running visualizations..." in text or "visualizations" in text.lower():
                for filename in pdb_filenames:
                    if filename in st.session_state.file_status:
                        if st.session_state.file_status[filename].get('minimization') == "✅":
                            # Only mark for visualization if minimization was successful
                            if 'visualization' not in st.session_state.file_status[filename] or st.session_state.file_status[filename]['visualization'] == '':
                                st.session_state.file_status[filename]['visualization'] = "In Queue"
                                print(f"  -> Marked {filename} as In Queue (visualization)")
            
            # Update the status display
            with status_container.container():
                # Show Minimization Phase
                st.markdown("### Minimization Phase")
                for filename in pdb_filenames:
                    if filename in st.session_state.file_status:
                        status = st.session_state.file_status[filename]
                        min_status = status.get('minimization', 'In Queue')
                        
                        if min_status == "processing":
                            st.markdown(f"**{filename}** ⚙️ Processing")
                        elif min_status == "✅":
                            st.markdown(f"**{filename}** ✅ Completed!")
                        else:
                            st.markdown(f"**{filename}** ⏳ In Queue")
                
                # Show Visualization Phase if any files have visualization status
                show_viz_section = False
                for filename in pdb_filenames:
                    if filename in st.session_state.file_status:
                        viz_status = st.session_state.file_status[filename].get('visualization', '')
                        if viz_status != '':
                            show_viz_section = True
                            break
                
                if show_viz_section:
                    st.markdown("---")
                    st.markdown("### Visualization Phase")
                    for filename in pdb_filenames:
                        if filename in st.session_state.file_status:
                            viz_status = st.session_state.file_status[filename].get('visualization', '')
                            
                            if viz_status == "processing":
                                st.markdown(f"**{filename}** ⚙️ Processing")
                            elif viz_status == "✅":
                                st.markdown(f"**{filename}** ✅ Completed!")
                            elif viz_status == "In Queue":
                                st.markdown(f"**{filename}** ⏳ In Queue")
        
        try:
            # Reset processing state
            st.session_state.processing_stopped = False
            st.session_state.processed_files = []
            st.session_state.visualization_images = []
            st.session_state.current_progress = []

            # Process all proteins directly in the specified output directory
            update_progress("Starting energy minimization...", 5)
            
            try:
                all_energy_data, protein_pairs = process_all_proteins(
                    pdb_files, 
                    output_path,  # Use user-specified output directory
                    selected_ff, 
                    selected_sm, 
                    iterations, 
                    threads, 
                    use_gpu,
                    None,  # gpu_platform_name
                    update_progress
                )
            except Exception as e:
                update_progress(f"Error in processing: {str(e)}", 100)
                st.session_state.processing_stopped = True
                raise
            
            # Save energy data
            update_progress("Saving energy data...", 90)
            energies_file = save_energy_data(all_energy_data, output_path)
            
            if run_visualizations and protein_pairs:
                update_progress("Running visualizations...", 92)

                # Import write_combined_rmsd_csv from utils
                from utils import write_combined_rmsd_csv

                # Pass image parameters to run_rmsd_analysis
                all_rmsd_data, successful_rmsd = run_rmsd_analysis(
                    protein_pairs, 
                    output_path,
                    update_progress,
                    image_width,
                    image_height,
                    image_dpi
                )
                
                if all_rmsd_data:
                    update_progress("Saving visualization data...", 98)
                    rmsd_csv, summary_csv = write_combined_rmsd_csv(all_rmsd_data, output_path)
            
            # Save metadata
            params = {
                'forcefield_info': selected_ff,
                'solvent_info': selected_sm,
                'iterations': iterations,
                'threads': threads,
                'run_rmsd': run_visualizations
            }
            metadata_file = save_metadata(output_path, params, all_energy_data)
            
            update_progress("Processing complete!", 100)

            # Mark processing as complete
            st.session_state.processing_stopped = True
            
            # Store results in session state
            st.session_state.results = {
                'energy_data': all_energy_data,
                'output_dir': str(output_path),
                'pdb_files': [str(f) for f in pdb_files]
            }
            
            # Display results
            st.markdown("## 📈 Results")
            
            # Success statistics
            successful_runs = [d for d in all_energy_data if d['initial_energy'] != 'ERROR']
            failed_runs = len(all_energy_data) - len(successful_runs)
            
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Total Files", len(all_energy_data))
            with col2:
                st.metric("Successful", len(successful_runs))
            with col3:
                st.metric("Failed", failed_runs)

            # Display visualization images if available
            if run_visualizations:
                viz_dir = output_path / "visualizations"
                if viz_dir.exists():
                    image_files = list(viz_dir.glob("*.png"))
                    if image_files:
                        st.markdown("### 📸 Visualization Images")

                        # Create two columns
                        col1, col2 = st.columns(2)

                        for idx, img_path in enumerate(image_files):
                            protein_name = img_path.stem  # filename without extension

                            if idx % 2 == 0:
                                with col1:
                                    st.image(
                                        str(img_path),
                                        caption=protein_name,
                                        width=500
                                    )
                            else:
                                with col2:
                                    st.image(
                                        str(img_path),
                                        caption=protein_name,
                                        width=500
                                    )

                    else:
                        st.info("No visualization images found.")

        except Exception as e:
            st.session_state.processing_stopped = True
            st.error(f"❌ Error during processing: {str(e)}")
            st.error(traceback.format_exc())
        finally:
            # Mark processing finished, but DO NOT clear state
            st.session_state.processing = False


if __name__ == "__main__":
    main()