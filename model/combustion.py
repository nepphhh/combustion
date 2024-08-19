import numpy as Math
import solver
from nist import elementalComposition
from nist import getElementsInSpecies as aij
from nist import getEnthalpyAtT as h
from nist import getEntropyAtT as s

# Solver parameters
precision = 1e-15

# System parameters
R = 8.31446261815324 		# J / K / mol
temperature = T0 = 3000. 	# K
pressure = P0 = 32.			# bar
composition = {				# in mols
	'H2': 3.17,
	'O2': 1, 
}

# Start tracking all the things we've got to track
elements, gaseous, condensed = set(), set(), set()

# Load up the elements automagically
for formula in composition:
	for element in elementalComposition[formula]:
		elements.add(element)
elements = list(elements)
elements.sort()

# Load up the species automagically
for formula in elementalComposition:
	if Math.all([element in elements for element in elementalComposition[formula]]):
		gaseous.add(formula)
gaseous, condensed = list(gaseous), list(condensed)
gaseous.sort()
condensed.sort()

# Get the full list of species
species = gaseous + condensed
I, J = len(elements), len(species)

# Initial element amounts
elementMols = {}
for i in elements:
	elementMols[i] = Math.sum([composition[j] * aij(i, j) for j in composition])

# Start tracking iteration variables
mols, lagrange_is= {}, {}
for j in species:
	mols[j] = 1
for i in elements:
	lagrange_is[i] = 0

# Iterate!
print("Beginning iteration...")
proceed, iteration, damps = True, 0, [1 for i in range(J)]
while proceed: 

	# Solve the matrix
	try:
		x = solver.solve(temperature, pressure, elements, gaseous, condensed, mols, elementMols, lagrange_is)
	except ValueError as e:
		print(e, 'Iteration halted at i =', iteration)
		break

	# Update variables
	for index in range(J):
		while mols[species[index]] + x[index] * damps[index] <= 0:
			damps[index] /= 2
		mols[species[index]] += x[index] * damps[index]
		if damps[index] != 1:
			damps[index] *= 2
	for index in range(I):
		lagrange_is[elements[index]] += x[index+J]

	# Keep going? 
	iteration += 1
	proceed = Math.any([Math.abs(delta) > precision for delta in x[:J]])

	# Print thingy
	if Math.mod(iteration, 1000)==0:
		print("...n="+str(iteration))

# A function (NOTE: returns in kJ NOT J!!!)
def gibbs(comp): 
	tot = 0
	n = 0
	for j in comp:
		if j in gaseous:
			n += comp[j]
	for j in comp:
		potiential = h(j, temperature) - temperature * s(j, temperature)
		if j in gaseous:
			potiential += R * temperature * Math.log(pressure * comp[j] / n) 
		tot += potiential * comp[j]
	return tot*1e-4

# Returns in J/K
def entropy(comp): 
	tot = 0
	n = 0
	for j in comp:
		if j in gaseous:
			n += comp[j]
	for j in comp:
		S = s(j, temperature)
		if j in gaseous:
			S -= R * Math.log(pressure * comp[j] / n) 
		tot += S * comp[j]
	return tot

# Clean up & sort
sortedSpecies = [{'name': j, 'amount': mols[j]} for j in species]
sortedSpecies.sort(key=lambda e: -e['amount'])

# Interface
print('\n- Iteration Complete (n='+str(iteration)+') -')
print('T:', temperature, 'K')
print('P:', pressure, 'bar')

print('\n- Reactant Composition -')
for j in composition:
	print(j+':', composition[j], 'mol')

totalMols = Math.sum([j['amount'] for j in sortedSpecies])
print('\n- Product Composition -')
for val in sortedSpecies:
	if val['amount'] > precision:
		print((val['name']+':').ljust(8), val['amount']/totalMols, 'mol fraction')
print('Total:', totalMols, 'mol')

print('\n- Delta G, Delta S -')
print(gibbs(mols) - gibbs(composition), 'kJ')
print(entropy(mols) - entropy(composition), 'J/K')


print('\n- Elemental Residuals -')
for i in elements:
	residual = Math.sum([mols[j] * aij(i, j) for j in species]) - elementMols[i]
	print(i+':', residual, 'mol')
