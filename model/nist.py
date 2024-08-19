import urllib.request
import csv, os

print("Obtaining NIST thermodynamic data...")

thermodynamicData = dict()

# URL References
NISTCodes = {

	'H':	'H-001',
	'NH':	'H-026', 
	'HNO':	'H-027',
	'HNO2 (Cis)':	'H-028',
	'HNO2 (Trans)':	'H-029',
	'HNO3':	'H-030',
	'OH':	'H-038',
	'HO2':	'H-043',
	'H2':	'H-050',
	'NH2':	'H-060', 
	'H2O':	'H-064',
	# 'H2O2':	'H-070',
	'NH3':	'H-083', 
	'N2H4': 'H-091', 

	'O':	'O-001',
	'O2':	'O-029',
	'O3':	'O-056',

	'CO':	'C-093',
	'CO2':	'C-095',
	'CH':	'C-043',
	'HCN':	'C-052',
	'HNCO': 'C-053', 
	'HCO':	'C-054',
	'CH2':	'C-057',
	'H2CO':	'C-061',
	'CH3':	'C-062',
	'CH4':	'C-067',
	'CN':	'C-080', 
	'CNO':	'C-087', 
	'CN2':	'C-088', 
	'C2':	'C-113',
	'C2H':	'C-124', 
	'C2H2': 'C-127', 
	'C2H4': 'C-128', 
	'C2H4O': 'C-129',
	'C2N': 	'C-133',
	'C2N2': 'C-134',
	'C3':	'C-138',

	'N':	'N-002', 
	'NO':	'N-005', 
	'NO2':	'N-007', 
	'NO3':	'N-009', 
	'N2':	'N-023', 
	'N2O':	'N-026', 
	'N3':	'N-034', 

}

def parseCompositionFromEquation(code):
	composition = {}
	lastchar = ''
	for char in code:
		if char == ' ':
			break
		if char == '-':
			composition['e-'] = 1
			continue
		if char == '+':
			composition['e-'] = -1
			continue
		if char.isalpha():
			composition[char] = 1
			lastchar = char
			continue
		if char.isdigit():
			composition[lastchar] = int(char)
			continue

	return composition

# This will be used in many other scripts
elementalComposition = {code: parseCompositionFromEquation(code) for code in NISTCodes}

# Load enthalpies from the web
for species in NISTCodes:

	code, file = NISTCodes[species], os.getcwd() + '/nist/'+species+'.csv'
	
	# Data structure for the data
	thermodynamicData[species] = dict()

	# Obtain from remote site
	if not os.path.isfile('filename.txt'):
		url = 'https://janaf.nist.gov/tables/' + code + '.txt'
		urllib.request.urlretrieve(url, file)

	# Open & parse
	with open(file) as csvFile:

		# Open file
		csvReader = csv.reader(csvFile, delimiter='\t')

		# Skip the first two rows
		rowskip = 2

		# Extract enthalphy information about each species
		for row in csvReader:

			# Skip the first two rows
			if rowskip: 
				rowskip -= 1
				continue

			# Reparse their infinity name
			if row[7] == 'INFINITE':
				row[7] = 'Infinity'
			
			# Store values
			try:
				thermodynamicData[species][float(row[0])] = {
					'Cp': float(row[1]), # specific heat at constant pressure
					'S': float(row[2]), # entropy (J/mol/K)
					'H': float(row[4])*1000, # enthalpy (J/mol)
					'Hf': float(row[5])*1000, # enthalpy of formation (J/mol)
				}
			except ValueError:
				if (row[0] == ''):
					continue
				else:
					print(row)
					exit('Critical problem parsing data for species ' + species + '.')


# What this package is all about: a method to retrieve interpolated thermo data
def getThermodynamicData(species, temperature): 
	
	# Load all temperatures
	temperatures = [t for t in thermodynamicData[species]]

	# Get the coldest temperature we have
	index = 0
	for i in range(len(temperatures)):
		if temperature < temperatures[i]:
			break
		index += 1

	# These are the temperatures to interpolate between
	floorTemp = temperatures[index-1]
	ceilTemp = temperatures[index]

	# Interpolation value
	interp = (temperature-floorTemp) / (ceilTemp-floorTemp)

	# Prepare to output
	output = {}

	# Get interpolated values
	for key in thermodynamicData[species][floorTemp]:
		floorVal = thermodynamicData[species][floorTemp][key]
		ceilVal = thermodynamicData[species][ceilTemp][key]
		output[key] = (1-interp) * floorVal + interp * ceilVal

	# Wrap up
	return output

def getElementsInSpecies(element, species):
	try:
		val = elementalComposition[species][element]
	except KeyError:
		val = 0
	return val

def getEntropyAtT(species, temperature):
	try:
		return getThermodynamicData(species, temperature)['S']
	except IndexError:
		exit("Error finding entropy of "+species+" at T="+str(temperature))

def getEnthalpyAtT(species, temperature):
	try:
		return getThermodynamicData(species, 298.15)['Hf'] + getThermodynamicData(species, temperature)['H']
	except IndexError:
		exit("Error finding enthalpy of "+species+" at T="+str(temperature))

def getSpecificHeatAtT(species, temperature):
	try:
		return getThermodynamicData(species, temperature)['Cp']
	except IndexError:
		exit("Error finding specific heat of "+species+" at T="+str(temperature))