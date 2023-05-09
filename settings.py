import json

# Load the JSON data from a file
with open('CELL_input.json') as file:
    data = json.load(file)

# Access the entries in the JSON data
packmol_entries = data['PACKMOL']
packmol_tolerance = packmol_entries['packmol_tolerance']
cell_radius = packmol_entries['cell_radius']
initial_packing_shape = packmol_entries['initial_packing_shape']
number_of_cells = data['number_of_cells']
beads = data['beads']

# Print the values of the accessed entries
print("Packmol Tolerance:", packmol_tolerance['value'])
print("Cell Radius:", cell_radius['value'])
print("Initial Packing Shape:", initial_packing_shape['value'])
print("Number of Cells:", number_of_cells['value'])

# Access individual beads
bead_A = beads['A']
bead_B = beads['B']
bead_N = beads['N']

print("Bead A Quantity:", bead_A['quantity'])
print("Bead A mass:", bead_A['mass'])
print("Bead B Quantity:", bead_B['quantity'])
print("Bead B mass:", bead_B['mass'])
print("Bead N Quantity:", bead_N['quantity'])
print("Bead N mass:", bead_N['mass'])

#Read JSON through a loop
with open('your_file.json') as file:
    data = json.load(file)

for category, subcategories in data.items():
    print(category)
    for subcategory, subcategory_data in subcategories.items():
        print(f" - {subcategory}")
        for key, value in subcategory_data.items():
            print(f"   {key}: {value}")
