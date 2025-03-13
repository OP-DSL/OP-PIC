import os
import re
import csv
from collections import defaultdict

def extract_route_values_to_csv(parent_folder, route_names, output_csv):
    # Regex pattern to extract G, M, D, and ARR values from filenames
    file_pattern = re.compile(r"log_G(\d+)_M(\d+)_D([\d.eE+-]+)_ARR(\d+)_R\d+\.log")

    # Dictionary to store minimum values grouped by (G, M, D, ARR)
    grouped_results = defaultdict(lambda: {route: float('inf') for route in route_names})

    # Walk through all subfolders and files
    for root, _, files in os.walk(parent_folder):
        for file_name in files:
            file_path = os.path.join(root, file_name)

            # Match the file pattern
            match = file_pattern.match(file_name)
            if not match:
                continue

            # Extract G, M, D, and ARR
            g_value = int(match.group(1))
            m_value = int(match.group(2))
            d_value = float(match.group(3))
            arr_value = int(match.group(4))

            # Read the log file
            with open(file_path, 'r') as file:
                lines = file.readlines()

            # Update the minimum values for each route
            for route in route_names:
                for line in lines:
                    if route in line:
                        parts = line.split()
                        if len(parts) > 4:
                            value = float(parts[4])
                            grouped_results[(g_value, m_value, d_value, arr_value)][route] = min(
                                grouped_results[(g_value, m_value, d_value, arr_value)][route],
                                value
                            )
                        break

    # Prepare results for CSV output
    sorted_results = sorted(grouped_results.items(), key=lambda x: (x[0][0], x[0][1], x[0][2], x[0][3]))

    # Write to CSV
    with open(output_csv, 'w', newline='') as csv_file:
        fieldnames = ["G", "M", "D", "ARR"] + route_names
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()

        for (g_value, m_value, d_value, arr_value), route_values in sorted_results:
            row = {
                "G": g_value,
                "M": m_value,
                "D": f"{d_value:.2e}",
                "ARR": arr_value,
            }
            for route in route_names:
                row[route] = (
                    route_values[route] if route_values[route] != float('inf') else None
                )
            writer.writerow(row)

# Input parameters
parent_folder = "./"  # Replace with the parent folder containing subfolders
route_names = ["MainLoop",
               "calculate_new_pos_vel_kernel", 
               "compute_electric_field_kernel", 
               "compute_node_charge_density_ke", 
               "compute_phi", 
               "deposit_charge_on_nodes_kernel", 
               "inject_ions_kernel", 
               "move_kernel",
               "Mv_Finalize",
               "Setup"]
output_csv = "output_sorted_routes.csv"

# Run the script
extract_route_values_to_csv(parent_folder, route_names, output_csv)
