import collections

#First we define the dictionary by name
#with no entries to start with
aerosol_dict = {}

#Now we can add a few entries
aerosol_dict['type'] = 'Sodium Chloride'
aerosol_dict['size [microns]'] = 10.56
aerosol_dict['Measured by'] = 'Owen'
aerosol_dict['values'] = [103,304,112,33]

# Cycle through each key:value pair and
# print to the screen
for key, value in aerosol_dict.items():
    print("Key  = ", key)
    print("Value = ", value)

# We can also check whether a key exists
if 'type' in aerosol_dict.keys():
  print("It is in there!")
else:
  print("Not there")

#We can also create 2D dictionaries which starts
#to become useful when embedding a structure to
#our key:value pairs. For example, in the following
#we manually create a 2D dictionary of aerosol type
#(by composition) and concentrations measured in
#a particular bin structure
aerosol_records = {'Sodium Chloride' : {'size [microns]':10.56, \
 'place':'coast', 'values':[103,304,112,33]}, \
 'Ammonium Sulphate' : {'size [microns]':0.345, \
  'place':'city', 'values':[450,1500,2003,579]}}
#We can then 'look up' and print the concentration array
#for ammonium sulphate as follows:
print(aerosol_records['Ammonium Sulphate']['values'])

# An alternative approach is to use the collections
#module. This has advantages over the previous approach
#if you want to avoid error messages with missing values
#or retain the order in which entries were added to
#the dictionary. For example the following takes an
#existing list and is able to append values to existing
#keys
aerosol_list = [('NaCl', 10.3), ('NH4NO3', 0.45), \
    ('NaCl', 23.5), ('NH4NO3', 0.34), ('NaCl', 12.1)]
aerosol_dict2 = collections.defaultdict(list)
for type, value in aerosol_list:
    aerosol_dict2[type].append(value)
print(aerosol_dict2.items())
