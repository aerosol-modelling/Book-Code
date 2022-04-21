# Import the relevant modules
import re
import pdb
import collections

# First read in the text file and print extracted equations to screen
filename= 'example_mechanism.txt'
text=open(filename,'rU')
with open (filename, "r") as myfile:
    data=myfile.read().replace('\n', '').replace('\t', '')

equation_numbers = re.findall(r"\{(.*?).\}",data)
print("Equation numbers = ",equation_numbers)
print("Total number of equations = ",len(equation_numbers))
eqn_list=re.findall(r"\}(.*?)\;",data)
print("Equation list = ",eqn_list)

# Initialise dictionaries that will store information on the loss and gain
# reactions for each compound. Also create dictionaries that store the
# stoichiometry and rate coefficient expressions as strings. Finally
# we create a dictionary that will allow us to map a compound name to
# an index that is used in the numerical arrays of our simulation.

rate_dict=collections.defaultdict(
          lambda: collections.defaultdict())
loss_dict=collections.defaultdict(
          lambda: collections.defaultdict())
gain_dict=collections.defaultdict(
          lambda: collections.defaultdict())
stoich_dict=collections.defaultdict(
          lambda: collections.defaultdict())
rate_dict_reactants=collections.defaultdict(
          lambda: collections.defaultdict())
species_dict=collections.defaultdict()
species_dict2array=collections.defaultdict()

#Create an integer that stores number of unique species
species_step=0

# In this loop we interrogate each line that has been extracted from our
# mechanism file and is now stored as a list
for equation_step in range(len(equation_numbers)):
    equation_full=eqn_list[equation_step]
    #split the line into reactants and products
    equation=equation_full.split(':',1)[0].split('=',1)
    # extract content to the left of the previous split [reactants]
    reactants=equation[0].split('+')
    #strip away all whitespace
    reactants= [x.strip(' ') for x in reactants]
    # extract content to the right of the previous split [products]
    products=equation[1].split('+')
    #strip away all whitespace
    products = [x.strip(' ') for x in products]
    #At the moment, we have not separated the reactant/product from
    #its stoichiometric value
    # Now extract the reaction rate expression
    rate_full=equation_full.split(':',1)[1]
    #strip away all whitespace
    rate_full=rate_full.strip()
    rate_dict[equation_step]="".join(rate_full.split())

    #used to identify reactants by number, for any given reaction
    reactant_step=0
    product_step=0

    # Loop through 'reactants' in this equation, extracting the
    # stoichiometry
    for reactant in reactants:
        reactant=reactant.split()[0]
        # - Extract stoichiometry and unique identifier
        try: #Extract a stoichiometric coefficient, if given
            temp=re.findall(r"[-+]?\d*\.\d+|\d+|\d+",reactant)
            #This extracts all numbers either side
            stoich=temp[0] #Select the first number extracted.
            # Check if this value is before the variable
            # If after, we ignore this. EG. '2NO2'
            # or just 'NO2'
            # If len(temp)==1 then we only have one number
            # and can proceed with the following
            if len(temp)==1:
                if reactant.index(stoich) == 0 :
                    reactant=reactant.split(stoich,1)[1]
                    stoich=float(stoich)
                else:
                    stoich=1.0
            elif len(temp)>1:
                # If this is the case, ensure the reactant
                # extraction is unique. If string is '2NO2'
                # the above procedure extracts 'NO'.
                # We need to ensure the reactant is 'NO2'.
                # We cut the value in temp[0] away from the
                # original string. We can attach the first
                # part with the second number. Thus
                if reactant.index(stoich) == 0 :
                    reactant=reactant.split(stoich,1)[1]+temp[1]
                    stoich=float(stoich)
                else:
                    stoich=1.0
        except:
            stoich=1.0
        # Store stoichiometry and species flags in dictionaries
        if reactant not in ['hv']:
            stoich_dict[equation_step][reactant_step]=stoich
            rate_dict_reactants[equation_step][reactant_step]=reactant
            # -- Update species dictionaries --
             #check to see if entry already exists
            if reactant not in species_dict.values():
                species_dict[species_step]=reactant
                species_dict2array[reactant]=species_step
                species_step+=1
            # -- Update loss dictionaries --
            if equation_step in loss_dict[reactant]:
                loss_dict[reactant][equation_step]+=stoich
            else:
                loss_dict[reactant][equation_step]=stoich

        reactant_step+=1
    # Now loop through the 'products' in this equation, extracting the
    # stoichiometry
    if len(products) > 0:
        for product in products:
            try:
                 #remove all tables, newlines, whitespace
                product=product.split()[0]
		        # - Extract stochiometry and unique product identifier
                try: #Now try to extract a stoichiometric coefficient if given
                     #This extracts all numbers either side
                    temp=re.findall(r"[-+]?\d*\.\d+|\d+|\d+",product)
                    #This selects the first number extracted, if any.
                    stoich=temp[0]
                    # Check if this value is before the variable
                    # If after, we ignore this. EG. '2NO2'
                    # or just 'NO2'
                    # If len(temp)==1 then we only have one number
                    # and can proceed with the following
                    if len(temp)==1:
                        if product.index(stoich) == 0 :
                            product=product.split(stoich,1)[1]
                            stoich=float(stoich)
                        else:
                            stoich=1.0
                    elif len(temp)>1:
                    # If this is the case, ensure the reactant
                    # extraction is unique. If string is '2NO2'
                    # the above procedure extracts 'NO'.
                    # We need to ensure the reactant is 'NO2'.
                    # We cut the value in temp[0] away from the
                    # original string. We can attach the first
                    # part with the second number. Thus
                        if product.index(stoich) == 0 :
                            product=product.split(stoich,1)[1]+temp[1]
                            stoich=float(stoich)
                        else:
                            stoich=1.0
                except:
                    stoich=1.0
                # - Store stoichiometry and species flags in dictionaries
                if product not in ['hv']:
                    # -- Update species dictionaries --
                    # check to see if entry already exists
                    if product not in species_dict.values():
                        species_dict[species_step]=product
                        species_dict2array[product]=species_step
                        species_step+=1
                    # -- Update loss dictionaries --
                    if equation_step in gain_dict[reactant]:
                        gain_dict[product][equation_step]+=stoich
                    else:
                        gain_dict[product][equation_step]=stoich

                product_step+=1
            except:
                pass

# Now print an example output from our new dictionaries for e.g. :
print("Gain dictionary entries for NAPINAOOH")
for equation_step in list(gain_dict['NAPINAOOH'].keys()):
    print("Equation number: ",equation_numbers[equation_step],
    " with stoichiometry", gain_dict['NAPINAOOH'][equation_step])
print("Gain dictionary entries for NAPINAO")
for equation_step in list(gain_dict['NAPINAO'].keys()):
    print("Equation number: ",equation_numbers[equation_step],
    " with stoichiometry", gain_dict['NAPINAO'][equation_step])
print("Loss dictionary entries for NAPINBO2")
for equation_step in list(loss_dict['NAPINBO2'].keys()):
    print("Equation number: ",equation_numbers[equation_step],
    " with stoichiometry", loss_dict['NAPINBO2'][equation_step])
print("Loss dictionary entries for A")
for equation_step in list(loss_dict['A'].keys()):
    print("Equation number: ",equation_numbers[equation_step],
    " with stoichiometry", loss_dict['A'][equation_step])

# Print the dependencies to screem to demonstrate potential sparsity of the
# mechanism
for species, species_step in species_dict2array.items():
    reactants_list_gain=[]
    reactants_list_loss=[]
    equations_list_gain=[]
    equations_list_loss=[]
    if species in loss_dict.keys():
        equations_list_loss=[num for num, stoich in loss_dict[species].items()]
        # Now extract all reactants in these equations as a list
        reactants_list_loss=[]
        for equation in equations_list_loss:
            reactants_list_loss=reactants_list_loss+[reactant for step, reactant in rate_dict_reactants[equation].items()]
        # Now turn this into a unique list
        reactants_list_loss=set(reactants_list_loss)
        reactants_list_loss=list(reactants_list_loss)

    if species in gain_dict.keys():
        equations_list_gain=[num for num, stoich in gain_dict[species].items()]
        # Now extract all reactants in these equations as a list
        reactants_list_gain=[]
        for equation in equations_list_gain:
            reactants_list_gain=reactants_list_gain+[reactant for step, reactant in rate_dict_reactants[equation].items()]
        # Now turn this into a unique list
        reactants_list_gain=set(reactants_list_gain)
        reactants_list_gain=list(reactants_list_gain)

    # Now combine these two lists and again merge into a unique one
    final_list=reactants_list_loss+reactants_list_gain
    final_list=set(final_list)
    final_list=list(final_list)

    print("Species dependent on ", species ," are:",final_list)
