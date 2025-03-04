import subprocess

def update_script(var1, var2):
    script_path = "#Feed in structures overhaul.py"
    
    with open(script_path, "r") as file:
        lines = file.readlines()
    
    # Modify the variables
    new_lines = []
    for line in lines:
        if line.startswith("cleav_prop ="):
            new_lines.append(f"cleav_prop = {var1}\n")
        elif line.startswith("cleav_prop_struct ="):
            new_lines.append(f"cleav_prop_struct = {var2}\n")
        else:
            new_lines.append(line)

    # Write the updated script back
    with open(script_path, "w") as file:
        file.writelines(new_lines)

cleav_prop_values = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]

cleav_prop_struct_values = [0.01, 0.009, 0.008, 0.007, 0.006, 0.005, 0.004, 0.003, 0.002, 0.001]

for var1, var2 in zip(cleav_prop_values, cleav_prop_struct_values):

        print(f"Updating script with cleav_prop={var1}, cleav_prop_struct={var2}" + "\n" )
        update_script(var1, var2)  # Modify the script
        subprocess.run(["python3", "#Feed in structures overhaul.py"]) 
         # Run the modified script, the above python3 will be modified to "py", "-3", if on windows
