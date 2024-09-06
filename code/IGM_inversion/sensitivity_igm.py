import os
import json
import subprocess
from shutil import copyfile, copy
import random

def generate_random_params():
    param1 = random.uniform(0, 2.4)
    param2 = random.uniform(0, 50)
    param3 = random.uniform(0, 0.5)
    # Add more parameters as needed
    return {"opti_velsurfobs_std": param1, "opti_thkobs_std": param2, "opti_usurfobs_std": param3}

def modify_params_file(params_file, params):
    with open(params_file, 'r') as f:
        data = json.load(f)
    data.update(params)
    with open(params_file, 'w') as f:
        json.dump(data, f, indent=4)

def run_igm_script():
    subprocess.run(['igm_run'])

def main():
    parameter_sets = []
    
    num_parameter_sets = 100  # Number of parameter sets to generate

    for _ in range(num_parameter_sets):
        params = generate_random_params()
        parameter_sets.append(params)

    input_files = ['thk.tif', 'icemask.tif', 'icemaskobs.tif',
                   'thkobs.tif', 'usurf.tif', 'usurfobs.tif',
                   'uvelsurfobs.tif', 'vvelsurfobs.tif','thkinit.tif']  # List of input files to copy


    #for idx, params in enumerate(parameter_sets[52:], start=53):
    for idx, params in enumerate(parameter_sets):
        output_folder = f'output_{idx}'
        print(output_folder)
        os.makedirs(output_folder, exist_ok=True)

        # Copy igm_run.py to output folder
        copyfile('igm_run.py', os.path.join(output_folder, 'igm_run.py'))
        copyfile('write_tif2.py', os.path.join(output_folder, 'write_tif2.py'))

        # Copy input files to output folder
        for input_file in input_files:
            copy(input_file, output_folder)

        # Copy params.json to output folder and modify it
        params_file = os.path.join(output_folder, 'params.json')
        copyfile('params.json', params_file)
        modify_params_file(params_file, params)

        # Run igm_run.py
        os.chdir(output_folder)  # Change directory to output folder
        run_igm_script()
        os.chdir('..')  # Change back to the original directory

if __name__ == "__main__":
    main()
