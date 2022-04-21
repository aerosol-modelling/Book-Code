from openbabel import pybel
pybel_object = pybel.readstring('smi','CC(=O)C1CC(C1(C)C)CC(=O)OC')

# Print attributes of our Molecule object
print("Molecular weight = ", pybel_object.molwt)
print("Molecular formula = ", pybel_object.formula)

# Print the
smarts_object = pybel.Smarts('[#6]')
print("SMARTS matches = ", smarts_object.findall(pybel_object))
print("Number of SMARTS matches = ", len(smarts_object.findall(pybel_object)))  
