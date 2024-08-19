import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy import log as ln
from numpy import cos, sin, pi, sqrt
import sys, math
from matplotlib.colors import hsv_to_rgb

# Functions
unit = lambda v: v / np.linalg.norm(v)
relative_angle = lambda v1, v2: np.arccos(np.clip(np.dot(unit(v1), unit(v2)), -1.0, 1.0))
relative_plane_angle = lambda p1, p2, v: np.abs(np.arccos(np.clip(np.dot(unit(np.cross(p1, p2)), unit(v)), -1.0, 1.0)) - (pi/2))

## DEFINE SYSTEM
# Function
f = lambda x, y, z: x+y+z
fmin, fmax = -8, 8

# Function derivative
df = lambda i, x, y, z: np.array([i, i, i])

# Constraints
c1 = lambda x, y, z: x*x - y*y - z 
c2 = lambda x, y, z: sqrt(x*x+y*y+z*z)-1

# Paramaterized constraints
c1_para = lambda i, j, k: np.array([j, k, j*j - k*k]) 
c2_para = lambda i, j, k: np.array([
	cos(j)*cos(k), sin(j)*cos(k), sin(k)]) 
map_view_para = lambda i, j, k: np.array([j, k, 2*i])
sphere_para = lambda i, j, k: np.array([
	cos(j)*cos(k), sin(j)*cos(k), sin(k)]) 

# Constraint gradients
dc1 = lambda i, x, y, z: np.array([2*x, -2*y, -i])
dc2 = lambda i, x, y, z: np.array([
	x/sqrt(x*x+y*y+z*z), 
	y/sqrt(x*x+y*y+z*z), 
	z/sqrt(x*x+y*y+z*z)
])

## DO CALCULATIONS
# Set projections
j1min, dj1, j1max = -2, 	0.02, 	2
k1min, dk1, k1max = -2, 	0.02, 	2
j2min, dj2, j2max = 0, 		0.02, 	2*pi
k2min, dk2, k2max = -pi/2, 	0.02, 	pi/2

# The grid of parameter values
j1, k1 = np.mgrid[slice(j1min, dj1+j1max, dj1), slice(k1min, dk1+k1max, dk1)]
j2, k2 = np.mgrid[slice(j2min, dj2+j2max, dj2), slice(k2min, dk2+k2max, dk2)]

# Get an identity set
i1 = np.ones(np.shape(j1))
i2 = np.ones(np.shape(j2))

# Get coordinates after projections
x1, y1, z1 = c1_para(i1, j1, k1)
x2, y2, z2 = c2_para(i1, j2, k2)

# Evaluate function in each projection
fn1 = f(x1, y1, z1)
fn2 = f(x2, y2, z2)

# Evalulate equality constraints (eqn = 0)
eqn1 = c2(x1, y1, z1)
eqn2 = c1(x2, y2, z2)

# Evaluate gradients, return as matrix of HSV color values
c1_grad1 = dc1(i1, x1, y1, z1)
c2_grad1 = dc2(i1, x1, y1, z1)
fn_grad1 = df(i1, x1, y1, z1)
deg1 = np.empty(np.shape(i1))
for n in range(np.shape(i1)[0]):
	for m in range(np.shape(i1)[1]):
		deg1[n,m] = relative_plane_angle(
			c1_grad1[:,n,m], c2_grad1[:,n,m], fn_grad1[:,n,m])
c1_grad2 = dc1(i2, x2, y2, z2)
c2_grad2 = dc2(i2, x2, y2, z2)
fn_grad2 = df(i2, x2, y2, z2)
deg2 = np.empty(np.shape(i2))
for n in range(np.shape(i2)[0]):
	for m in range(np.shape(i2)[1]):
		deg2[n,m] = relative_plane_angle(
			c1_grad2[:,n,m], c2_grad2[:,n,m], fn_grad2[:,n,m])

## PLOT FUNCTION
# The colorplot we want
cmap0, cmap1 = 'RdBu', 'twilight'

# Set up subplots
fig, ((ax_fn1, ax_deg1), (ax_fn2, ax_deg2)) = plt.subplots(nrows=2, ncols=2)

# Plot function as colormesh with contours in the first projection
fn_plot1 = ax_fn1.pcolormesh(j1, k1, fn1, 
	cmap=cmap0, shading='auto', vmin=fmin, vmax=fmax)
ax_fn1.contour(j1, k1, fn1, 
	cmap=cmap0, vmin=fmin, vmax=fmax)
fig.colorbar(fn_plot1, ax=ax_fn1)
# Plot the equality equation
ax_fn1.contour(j1, k1, eqn1, levels=[0], cmap=cmap0)
# Add labels
ax_fn1.set_xlabel('j1')
ax_fn1.set_ylabel('k1')
ax_fn1.set_aspect('equal')
ax_fn1.xaxis.set_ticks(np.arange(j1min, j1max, 1))
ax_fn1.yaxis.set_ticks(np.arange(k1min, k1max, 1))

# Plot the gradient angle difference in the first projection
deg_plot1 = ax_deg1.pcolormesh(j1, k1, deg1, 
	cmap=cmap1, shading='auto', vmin=0, vmax=pi/2)
ax_deg1.contour(j1, k1, deg1, 
	cmap=cmap1, vmin=0, vmax=pi/2)
fig.colorbar(deg_plot1, ax=ax_deg1)
# Plot the equality equation
ax_deg1.contour(j1, k1, eqn1, levels=[0], cmap=cmap1)
# Add labels
ax_deg1.set_xlabel('j1')
ax_deg1.set_ylabel('k1')
ax_deg1.set_aspect('equal')
ax_deg1.xaxis.set_ticks(np.arange(j1min, j1max, 1))
ax_deg1.yaxis.set_ticks(np.arange(k1min, k1max, 1))

# Plot function as colormesh with contours in the second projection
fn_plot2 = ax_fn2.pcolormesh(j2, k2, fn2, 
	cmap=cmap0, shading='auto', vmin=fmin, vmax=fmax)
ax_fn2.contour(j2, k2, fn2, 
	cmap=cmap0, vmin=fmin, vmax=fmax)
fig.colorbar(fn_plot2, ax=ax_fn2)
# Plot the equality equation
ax_fn2.contour(j2, k2, eqn2, levels=[0], cmap=cmap0)
# Add labels
ax_fn2.set_xlabel('j2')
ax_fn2.set_ylabel('k2')
ax_fn2.set_aspect('equal')
ax_fn2.xaxis.set_ticks(np.arange(j2min, j2max, 1))
ax_fn2.yaxis.set_ticks(np.arange(k2min, k2max, 1))

# Plot the gradient angle difference in the second projection
deg_plot2 = ax_deg2.pcolormesh(j2, k2, deg2, 
	cmap=cmap1, shading='auto', vmin=0, vmax=pi/2)
ax_deg2.contour(j2, k2, deg2, 
	cmap=cmap1, vmin=0, vmax=pi/2)
fig.colorbar(deg_plot2, ax=ax_deg2)
# Plot the equality equation
ax_deg2.contour(j2, k2, eqn2, levels=[0], cmap=cmap1)
# Add labels
ax_deg2.set_xlabel('j2')
ax_deg2.set_ylabel('k2')
ax_deg2.set_aspect('equal')
ax_deg2.xaxis.set_ticks(np.arange(j2min, j2max, 1))
ax_deg2.yaxis.set_ticks(np.arange(k2min, k2max, 1))

# Set figure layout
fig.tight_layout(pad=0.1)

## LAGRANGE SPACE FUNCTION
# Choose point where function is to be evaluated in lagrange space
lag_x, lag_y, lag_z = 0.788095, 0.442572, 0.426423

# Evaluate partials at function position
dfdx, dfdy, dfdz = df(1, lag_x, lag_y, lag_z)
dc1dx, dc1dy, dc1dz = dc1(1, lag_x, lag_y, lag_z)
dc2dx, dc2dy, dc2dz = dc2(1, lag_x, lag_y, lag_z)

# Make equations for evaluation
lagrange_fn_x = lambda l1, l2: l1*dc1dx + l2*dc2dx - dfdx
lagrange_fn_y = lambda l1, l2: l1*dc1dy + l2*dc2dy - dfdy
lagrange_fn_z = lambda l1, l2: l1*dc1dz + l2*dc2dz - dfdz

# Get lagrange dimensions
l1min, dl1, l1max = -4, 0.2, 4
l2min, dl2, l2max = -4, 0.2, 4

# The grid of parameter values
l1, l2 = np.mgrid[slice(l1min, dl1+l1max, dj1), slice(l2min, dl2+l2max, dl2)]

# Evaluate functions
lagrange_x = lagrange_fn_x(l1, l2)
lagrange_y = lagrange_fn_y(l1, l2)
lagrange_z = lagrange_fn_z(l1, l2)

## PLOT THE LAGRANGE SPACE FIGURE
# Set up subplots
lagrange_fig, (lagrange_ax_x, lagrange_ax_y, lagrange_ax_z) = plt.subplots(ncols=3)

# Set up color maps & dimensions
lagrange_cmap = 'PiYG'
lagrange_min, lagrange_max = -10, 10

# Plot function as colormesh with contours in the first projection
lagrange_plot_x = lagrange_ax_x.pcolormesh(l1, l2, lagrange_x, 
	cmap=lagrange_cmap, shading='auto', vmin=lagrange_min, vmax=lagrange_max)
lagrange_ax_x.contour(l1, l2, lagrange_x, 
	cmap=lagrange_cmap, vmin=lagrange_min, vmax=lagrange_max)
lagrange_fig.colorbar(lagrange_plot_x, ax=lagrange_ax_x)
# Plot the equality equation
lagrange_ax_x.contour(l1, l2, lagrange_x, levels=[0], cmap=lagrange_cmap)
# Add labels
lagrange_ax_x.set_xlabel('l1')
lagrange_ax_x.set_ylabel('l2')
lagrange_ax_x.set_aspect('equal')

# Plot function as colormesh with contours in the first projection
lagrange_plot_y = lagrange_ax_y.pcolormesh(l1, l2, lagrange_y, 
	cmap=lagrange_cmap, shading='auto', vmin=lagrange_min, vmax=lagrange_max)
lagrange_ax_y.contour(l1, l2, lagrange_y, 
	cmap=lagrange_cmap, vmin=lagrange_min, vmax=lagrange_max)
lagrange_fig.colorbar(lagrange_plot_y, ax=lagrange_ax_y)
# Plot the equality equation
lagrange_ax_y.contour(l1, l2, lagrange_y, levels=[0], cmap=lagrange_cmap)
# Add labels
lagrange_ax_y.set_xlabel('l1')
lagrange_ax_y.set_ylabel('l2')
lagrange_ax_y.set_aspect('equal')

# Plot function as colormesh with contours in the first projection
lagrange_plot_z = lagrange_ax_z.pcolormesh(l1, l2, lagrange_z, 
	cmap=lagrange_cmap, shading='auto', vmin=lagrange_min, vmax=lagrange_max)
lagrange_ax_z.contour(l1, l2, lagrange_z, 
	cmap=lagrange_cmap, vmin=lagrange_min, vmax=lagrange_max)
lagrange_fig.colorbar(lagrange_plot_z, ax=lagrange_ax_z)
# Plot the equality equation
lagrange_ax_z.contour(l1, l2, lagrange_z, levels=[0], cmap=lagrange_cmap)
# Add labels
lagrange_ax_z.set_xlabel('l1')
lagrange_ax_z.set_ylabel('l2')
lagrange_ax_z.set_aspect('equal')

# Set figure layout
lagrange_fig.tight_layout(pad=0.1)

# Show thing
plt.show()
