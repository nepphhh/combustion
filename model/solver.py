import numpy as Math
from numpy import log as ln
from nist import getElementsInSpecies as aij
from nist import getEnthalpyAtT as h
from nist import getEntropyAtT as s

R = 8.31446261815324 

# Gaussian scaled partial pivot
def gaussian(A, b):
	
  n = len(b)
  
  # Keeps track of row order
  L = [i for i in range(n)]
  # Used for scaling the rows
  S = [0.0 for i in range(n)]
  # Solving for x vect
  x = [0.0 for i in range(n)]

  # Find scaling values
  for i in range(n-1):
    for j in range(n-1):
      S[i] = max(S[i], abs(A[i][j]))

  # Forward Elimination
  for k in range(n-1):
    R = 0.0
    for i in range(k, n-1):
      temp = abs(A[L[i]][k] / S[L[i]])
      if temp > R:
        R = temp
        index = i

    # Do the swap
    temp = L[index]
    L[index] = L[k]
    L[k] = temp

    # Zero out below the diagonal
    for i in range (k+1, n):
      xmult = A[L[i]][k]/A[L[k]][k]
      for j in range(k+1, n):
        A[L[i]][j] = A[L[i]][j] - xmult * A[L[k]][j]
      b[L[i]] = b[L[i]] - xmult * b[L[k]]

  # Back substitution
  x[n-1] = b[L[n-1]] / A[L[n-1]][n-1]

  for k in range(n-1, -1, -1):
    sum = b[L[k]]
    for j in range(k+1, n):
      sum = sum - A[L[k]][j] * x[j]
    x[k] = sum / A[L[k]][k]

  return x

# The iterator function
def solve (T, P, elements, gaseous, condensed, moles, b_0s, lambda_is):

	species = gaseous+condensed

	# Make some lil functions
	l = lambda i: lambda_is[i]
	b = lambda i: b_0s[i]
	n = lambda j: moles[j]
	N = Math.sum([n(j) for j in gaseous])

	if Math.min([n(j) for j in species]) <= 0:
		raise ValueError('Negative or zero mols of something.')

	# Build the iteration matrix
	A = (
		([
			[aij(i_, j) for j in condensed] +
			[aij(i_, j) for j in gaseous] +
			[0 for i in elements]
			for i_ in elements
		]) + ([
			[0 for j in condensed] + 
			[R*T/n(j) - R*T/N if j_==j else -R*T/N for j in gaseous] +
			[-aij(i, j_) for i in elements]
			for j_ in gaseous
		]) + ([
			[0 for j in condensed] + 
			[0 for j in gaseous] + 
			[-aij(i, j_) for i in elements]
			for j_ in condensed
		])
	)

	# Set up the right-hand-side
	b = (
		[b(i) - Math.sum([aij(i, j)*n(j) for j in species]) for i in elements] +
		[Math.sum([l(i)*aij(i, j) for i in elements]) - h(j, T) + T*s(j, T) - R*T*ln(P*n(j)/N) for j in gaseous] +
		[Math.sum([l(i)*aij(i, j) for i in elements]) - h(j, T) + T*s(j, T)for j in condensed]
	)

	# Solve
	return gaussian(A, b)
