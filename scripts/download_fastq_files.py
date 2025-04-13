import os
import pandas as pd
import subprocess
from pathlib import Path
from Bio import Entrez
import xml.etree.ElementTree as ET
import shutil
from profilehooks import profile

# Set your email here
Entrez.email = "spenser.johnson@yale.edu"

def fetch_srr_and_library_name(srx_id):
    handle = Entrez.efetch(db="sra", id=srx_id, rettype="xml")
    records = handle.read()
    root = ET.fromstring(records)
    srr_id_list = [run.attrib["accession"] for run in root.findall(".//RUN")]
    srr_id = srr_id_list[0] if len(srr_id_list) == 1 else srr_id_list
    library_name = root.find(".//LIBRARY_DESCRIPTOR/LIBRARY_NAME").text
    return srr_id, library_name
    
def download_srr(srr_id, library_name):
    os.system(f"prefetch {srr_id}")
    os.rename(srr_id, f"{library_name}_{srr_id}")

@profile(stdout=False, filename = r'C:\\Users\\spejo\\Documents\\1_CRISPR_analysis_test_input\\FASTQ\\download_fastq_files_baseline.prof')
def download_fastq_files(input_dir:str, input_file:str, overwrite:bool = False):
    # Read the input file
    df = pd.read_csv(input_file)
    
    # Assign file path to sra toolkit
    sra_toolkit_path = str(Path(r"C:\Users\spejo\Documents\Scripts\sratoolkit.3.1.1-win64\bin"))
    
    # Initialize variables
    output_data = {}

    for row in df.itertuples():
        # Create the directory path if it doesn't exist already
        dir_path = os.path.join(input_dir, row.First_Last, row.Screen_ID)

        if not Path(dir_path).exists():
            os.makedirs(dir_path)

        # Initialize treatment and control samples - will be used to make contrast table
        t_samples = ','.join([str(x).strip() for x in [row.treatment1, row.treatment2, row.treatment3] if pd.notnull(x)])
        c_samples = ','.join([str(x).strip() for x in [row.control1, row.control2, row.control3] if pd.notnull(x)])
        
        srx_ids = [str(x).strip() for x in [row.control1, row.control2, row.control3, row.treatment1, row.treatment2, row.treatment3] if pd.notnull(x)]

        for srx_id in srx_ids:
            srr_id, library_name = fetch_srr_and_library_name(srx_id)
            print(f"Found {srr_id} for {srx_id}")

            srr_id_exists = any(srr_id in str(file) for file in Path(dir_path).glob('*.fastq'))

            if srr_id_exists and not overwrite:
                pass
            else:
                prefetch_cmd = ["prefetch.exe",
                        srr_id,
                        "-O",
                        str(Path(dir_path))]
            
                prefetch_cmd = " ".join(prefetch_cmd)

                fasterq_dump_cmd = ["fasterq-dump.exe", 
                                    srr_id,
                                    "--split-files",
                                    "-O",
                                    str(Path(dir_path))]
                
                fasterq_dump_cmd = " ".join(fasterq_dump_cmd)
                
                os.chdir(sra_toolkit_path)
                print(f"Working directory set to: {sra_toolkit_path}")
                subprocess.run(prefetch_cmd, shell=True)
                print(f"Prefetched {srr_id}")
                
                os.chdir(sra_toolkit_path)
                subprocess.run(fasterq_dump_cmd, shell=True)
                print(f'Executed {fasterq_dump_cmd}')
                print(f"Converted {srr_id} .sra file into .fastq")

                # Remove the SRR directory
                shutil.rmtree(Path(dir_path) / Path(srr_id))
                print(f"Deleted directory {Path(dir_path) / Path(srr_id)}")

            new_name = f"{library_name}_{srx_id}_{srr_id}.fastq"
            new_path = Path(dir_path) / Path(new_name)

            if not Path(new_path).exists():
                os.rename(Path(dir_path) / f"{srr_id}_1.fastq", new_path)
                print(f"Renamed {Path(dir_path) / Path(srr_id)}_1.fastq to {new_path}")
                
            else:
                print(f"{new_path} already exists.")

            new_name_cnttbl = new_name.replace('.fastq', '')
            if srx_id in t_samples or srx_id in c_samples:
                t_samples = t_samples.replace(srx_id, new_name_cnttbl)
                c_samples = c_samples.replace(srx_id, new_name_cnttbl)

        # Collect output data for each Screen_ID
        if row.Screen_ID not in output_data:
            output_data[row.Screen_ID] = []
        output_data[row.Screen_ID].append([row.contrast, t_samples, c_samples])

    # Write the output data to files in their respective directories
    for screen_id, data in output_data.items():
        output_df = pd.DataFrame(data, columns=['contrast', 'treatment', 'control'])
        screen_dir = os.path.join(input_dir, row.First_Last, screen_id)
        output_df.to_csv(os.path.join(screen_dir, f"{screen_id}_cnttbl.txt"), sep='\t', index=False)

    print("All done!")

# Example usage
input_dir = r"C:\Users\spejo\Documents\1_CRISPR_analysis_test_input\FASTQ"
input_file = r"C:\Users\spejo\Documents\1_CRISPR_analysis_test_input\FASTQ\crispr_screen_bioprojects.csv"
input_path = Path(input_dir)

download_fastq_files(input_dir, input_file, overwrite=False)