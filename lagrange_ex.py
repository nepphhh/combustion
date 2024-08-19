import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy import log as ln
from numpy import cos, sin, arctan, sqrt, pi
import sys

# Define function
f = lambda x, y, z: cos(z)*sin(20*arctan(y/x)) + sin(z)*cos(5*ln(x*y)) # old, way too complex function
fmin, fmax = -1.5, 1.5
# f = lambda x, y, z: sin(x)+sin(y)+sin(z)
# fmin, fmax = -3, 3

# Define projection
proj = lambda t1, t2, k: (t1, t2, k) # The vertical projection
# proj = lambda t1, t2, k: ((t1-pi)**2 + (t2-pi)**2, t1, t2) # This is the paraboloid projection
# proj = lambda t1, t2, k: (pi + 1.7*cos(t2), t1, pi + 1.7*sin(t2)) # This is the cylindrical projection

# Define equations solved for
eq = lambda x, y, z: sqrt((x-pi)**2 + (z-pi)**2) - 1.7 # This is the intersection with the cylindrical constraint
# 	eq = lambda x, y, z: x - (y-pi)**2 - (z-pi)**2 # This is the intersection with the paraboloid constraint
# eq = lambda x, y, z:

# Generate a full grid
# dv, vmin, vmax = 1, 1, 3 # 0.01, 0.05, 3
# v = slice(vmin, dv+vmax, dv)
# x, y, z = np.mgrid[v, v, v]

# Set up however many projections we want
tmin, dt, tmax = 0, 0.01, 2*pi
T1, T2, Fn, Eqn, n = [], [], [], [], 25
for i in range(n):

	# The grid of parameter values
	t1, t2 = np.mgrid[slice(tmin, dt+tmax, dt), slice(tmin, dt+tmax, dt)]

	# An extra parameter, if needed
	k = i*0.2

	# Get coordinates after projection
	x, y, z = proj(t1, t2, k)

	# Evaluate function & constraint
	fn = f(x, y, z)
	eqn = eq(x, y, z)

	# Save values
	T1.append(t1)
	T2.append(t2)
	Fn.append(fn)
	Eqn.append(eqn)

# Make us a figure
def drawFigure(fig, i):
	ax = fig.add_subplot(111)
	cmap = plt.get_cmap('PiYG')
	colormesh = ax.pcolormesh(T1[i], T2[i], Fn[i], cmap=cmap, shading='auto', vmin=fmin, vmax=fmax)
	fn_contours = ax.contour(T1[i], T2[i], Fn[i], cmap=cmap, vmin=fmin, vmax=fmax)
	# eqn_contours = ax.contour(T1[i], T2[i], Eqn[i], levels=[0])
	fig.colorbar(colormesh)
	# ax.clabel(fn_contours, inline=1, fontsize=8)	
	plt.xlabel("t_1")
	plt.ylabel("t_2")
	plt.title("Projection " + str(i))
	fig.tight_layout(pad=1)

# Set up plot
i = 0
fig = plt.figure()
drawFigure(fig, i)

# What happens on key press
def onKeypressEvent(fig, key):
	sys.stdout.flush()
	global i
	if key == ',' and i > 0:
		i -= 1
	elif key == '.' and i < n-1:
		i += 1
	else:
		return 0
	fig.clear()
	drawFigure(fig, i)
	plt.draw()

# Start the key listener
fig.canvas.mpl_connect('key_press_event', lambda event: onKeypressEvent(fig, event.key))

# Show thing
plt.show()
