#%%
import pandas as pd
import numpy as np
import tkinter as tk
from tkinter import filedialog as fd
from io import StringIO
import re

#%%

print("This program reads bond lengths and angles from .cif files and saves them as a .txt to copy from.")
print("For this to work, the atoms need to be named [element][number][letter for disorders].")
print("Please check the output files, if all required data is included.")
input("Press Enter to continue...")
root = tk.Tk()
filename = fd.askopenfilename()
file_list = filename.split("/")
path = "/".join(file_list[:-1]) +"/"
root.destroy()

#%%

def extract_lines(start, lines):
    result = []
    for line in lines[start:]:
        if len(line) <= 3:
            break
        else:
            result.append(line)
    return result

file = open(filename, "r")
lines = file.readlines()
dist = []
angles = [] 

#############################################################################

dist_start = (lines.index("  _geom_bond_distance\n") + 3)
dist = extract_lines(dist_start, lines)

dist= "".join(dist)
dist_table = pd.read_csv(StringIO(dist), delimiter = " ", header = None)
dist_table = dist_table.drop(columns = [0, 4, 5])
dist_table.columns = ["Atom A", "Atom B", "Distance"]

print(dist_table)

#############################################################################

angle_start = (lines.index("  _geom_angle\n") + 4)
angles = extract_lines(angle_start, lines)

angles = "".join(angles)
angles_table = pd.read_csv(StringIO(angles), delimiter = " ", header = None)
angles_table = angles_table.drop(columns = [0, 5, 6, 7])
angles_table.columns = ["Atom A", "Atom B", "Atom C", "Angle"]

print(angles_table)

# %%

def is_in_elements(value):
    clean_value = []
    for character in value:
        if character.isnumeric():
            break
        clean_value.append(character)
        
    clean_value = "".join(clean_value)
    for element in elements:
        if clean_value == element:
            return True
    return False

def clean_element(value):
    clean_value = []
    for character in value:
        if character.isnumeric():
            break
        clean_value.append(character)
        
    clean_value = "".join(clean_value)
    return clean_value

def disorder(value):
    if value[-1].isnumeric():
        return None
    else: 
        return value[-1]

def same_disorder(rows):
    disorders = [disorder(r) for r in rows]
    none_count = sum(d == None for d in disorders)

    if none_count >= len(rows) - 1:
        return True
    
    x_num = None
    for x in disorders:
        if x != None:
            x_num = x
            break
    for x in disorders:
        if x != x_num and x != None:
            return False
    return True

elements = input("Bond of interest (e.g. C-H): ").split("-")

dist_roi = []

while "".join(elements) != "stop":
    atom_a = elements[0]
    atom_b = elements[1]

    for _,row in dist_table.iterrows():
        d0 = disorder(row[0])
        d1 = disorder(row[1])

        if not same_disorder([row[0], row[1]]):
            continue

        c0 = clean_element(row[0])
        c1 = clean_element(row[1])

        if (
        (atom_a == c0 and atom_b == c1)
        or (atom_b == c0 and atom_a == c1)
        ):
            dist_roi.append(row)

    print(f"Bond {'-'.join(elements)} added. Type another bond or type stop if done.")
    elements = input("Bonds of interest (e.g. C-H): ").split("-")
                   
dist_roi_table = pd.DataFrame(dist_roi)
dist_roi_string = dist_roi_table.to_string(header = False, index = False, index_names = False).split("\n")
dist_roi_string2 = []
for x in dist_roi_string:
    x2 = " ".join(x.split())
    dist_roi_string2.append(x2)
dist_roi_string = ", ".join(dist_roi_string2)

print(dist_roi_string)
print("Bonds done. Continue with angles.")

#######################################################################################

elements = input("Angles of interest (e.g. H-C-H): ").split("-")

angles_roi = []

while "".join(elements) != "stop":
    atom_a = elements[0]
    atom_b = elements[1]
    atom_c = elements[2]

    for _,row in angles_table.iterrows():
        if same_disorder([row[0], row[1], row[2]]):
            c0 = clean_element(row[0])
            c1 = clean_element(row[1])
            c2 = clean_element(row[2])
            if ((atom_a == c0 and atom_b == c1 and atom_c == c2)
            or (atom_a == c2 and atom_b == c1 and atom_c == c0)):
                angles_roi.append(row)

    print(f"Angle {'-'.join(elements)} added. Type another angle or type stop if done.")
    elements = input("Angles of interest (e.g. H-C-H): ").split("-")
                   
angles_roi_table = pd.DataFrame(angles_roi)
angles_roi_string = angles_roi_table.to_string(header = False, index = False, index_names = False).split("\n")
angles_roi_string2 = []
for x in angles_roi_string:
    x2 = " ".join(x.split())
    angles_roi_string2.append(x2)
angles_roi_string = ", ".join(angles_roi_string2)

print(angles_roi_string)

output_string = dist_roi_string + angles_roi_string

f = open("XRD Reader Output.txt", "w")
f.write(output_string)
file.close()
print("Output file generated.")
input("Press Enter to continue...")

# %%