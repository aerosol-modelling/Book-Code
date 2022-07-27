#------------------------------------------------------------------------------------
#
#  Calculate saturation vapour pressure of acetic acid and propanoic acid
#
#
#------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------
#
# 
# Import required packages
# 
#
#------------------------------------------------------------------------------------

import re
import pdb 
import pybel
import collections
import pkg_resources as pkg
from math import log
import glob


#------------------------------------------------------------------------------------
# 
# The following functions were created for UManSysProp and are given for the purpose
# of the exercise.
# 
# _read_data and _read_smarts from UManSysProp __init__.py
# matches and aggregate_matches from UManSysProp groups.py
#
#------------------------------------------------------------------------------------



def smarts(s):
    if not isinstance(s, bytes):
        s = s.encode('ascii')
    try:
        return pybel.Smarts(str(s,'utf-8'))
    except IOError as e:
        # Convert pybel's IOError (?!) into a ValueError
        raise ValueError(str(e))

_parse_re = re.compile(r'^\s*((?P<data>([^#]|\S#)*)(\s+#.*)?|#.*)$')
def _read_data(
        filename, key_conv=str, value_conv=float, key_col=0, value_col=1):
    result = {}
    cols = None
    for count, line in enumerate(
            pkg.resource_stream(__name__, filename), start=1):
        data = _parse_re.match(line.decode('utf-8')).group('data')
        if data:
            data = data.split()
            try:
                if cols is None:
                    cols = len(data)
                elif len(data) != cols:
                    raise ValueError(
                            'Unexpected number of values (expected %d)' % cols)
                key = key_conv(data[key_col])
                value = value_conv(data[value_col])
                if key in result:
                    raise ValueError(
                            'Duplicate definition for group %s' % key)
                result[key] = value
            except (IndexError, ValueError) as e:
                e.args += ('on line %d of %s' % (count, filename),)
                raise
    return result

def _read_smarts(filename):
    return _read_data(filename, key_conv=int, value_conv=smarts)

def matches(patterns, compound):
    """
    Returns a mapping of group identifier to number of matches.

    The *patterns* parameter specifies a dictionary mapping a group number to
    OpenBabel Smarts objects.

    The *compound* parameter specifies an OpenBabel Molecule object (presumably
    derived from a SMILES string). The function calculates the number of
    matches of the Molecule with each Smarts string, returning the result as a
    mapping.

    :param patterns:
        The mapping of group numbers to Smarts objects

    :param compound:
        A Molecule objects to match against the SMARTS strings in *filename*

    :returns:
        A mapping of group number to match counts
    """
    return {
        group: len(smarts.findall(compound))
        for group, smarts in patterns.items()
        }

def aggregate_matches(matches, coefficients, groups=None):
    """
    Calculates the biased sum of a set of matches.

    The *matches* parameter specifies a mapping of groups to values to be
    summed.  The *coefficients* parameter specifies mapping of groups to
    coefficient values.

    The optional *groups* parameter specifies the subset of the keys of
    *matches* which are to be summed. If this is not specified, it defaults to
    all keys of *matches*.

    For each group in *groups*, the corresponding values from *matches* and
    *coefficients* will be multiplied. The result is the sum of all these
    multiplications.

    :param matches:
        A mapping of groups to values to be summed

    :param coefficients:
        A mapping of groups to coefficients

    :param groups:
        An optional subset of the keys of *matches* indicating the values to
        be summed

    :returns:
        A floating point value which is the sum of the multiplication of each
        value in *matches* against the corresponding value in *coefficients*
    """
    if groups is None:
        groups = matches.keys()
    return sum(
        matches[group] * coefficients[group]
        for group in groups
        )

#------------------------------------------------------------------------------------
# define a variable that reads in the SMARTS
#
# YOUR CODE HERE:






#------------------------------------------------------------------------------------
# define function to search for groups
#
# YOUR CODE HERE:





#------------------------------------------------------------------------------------
# define variables for each parameter used in SIMPOL
#
# YOUR CODE HERE:





#------------------------------------------------------------------------------------
# define function to caculate vapour pressure from identified groups
#
# YOUR CODE HERE:





#------------------------------------------------------------------------------------

# The following code block is taken from the example folder of UManSysProp and can 
# be used to read the SMILES of multiple compounds using a single input.

#-------------------------------------------------------------------------------------

# Read in the property data, seperating only .prop files from the rest.
# I have defined the .prop files myself, just for ease of use for this project.
# It is just a textfile with compound name/reference as first column, SMILES as second.

onlyfiles = [f for f in glob.glob('*.prop')]

# extract the data from each file and start generating a dictionary of Pybel objects

step=0
filenames=[]
Compound_reference=[]
smiles_array=[]
property_array=[]
Pybel_object_dict=dict()

for filename in onlyfiles: # If you have more than one file, for whatever reason

   SMILES_flag=0
   filenames.append(filename[:])
   text=open(filename[:],'rU')
   
   for line in text:
       input = line.split()
       # Keep a list of the information
       Compound_reference.append(input[0])
       smiles_array.append(input[1])
       # Now create Pybel objects which are used in all property predictive techniques
       Pybel_object=pybel.readstring('smi',input[1])
       Pybel_object_dict[input[1]]=Pybel_object

vapour_pressure_dict=collections.defaultdict(lambda: collections.defaultdict())

#--------------------------------------------------------------------------------------
# Define a temperature to calculate vapour pressure of a compound at that temperature

temperature = 298.15

#--------------------------------------------------------------------------------------
#
# calculate vapour pressure at stated temperature [log10(atm)]
#
#--------------------------------------------------------------------------------------

for smiles in smiles_array:
    vapour_pressure_dict[smiles]['SIMPOL']=simpol_vapour_pressure(Pybel_object_dict[smiles], temperature)
    print(vapour_pressure_dict[smiles]['SIMPOL'])