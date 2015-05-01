"""
Useful things for plotting 1D GeoClaw results.
"""
import setrun
rundata=setrun.setrun()
drytol_default = rundata.geo_data.dry_tolerance

def topo(current_data):
   """
   topo is assumed to be aux[0,:]
   """
   aux = current_data.aux
   topo = aux[0,:]
   return topo

def land(current_data):
   """
   Return a masked array containing the surface elevation only in dry cells.
   """
   from numpy import ma
   drytol = getattr(current_data.user, 'dry_tolerance', drytol_default)
   q = current_data.q
   aux = current_data.aux
   h = q[0,:]
   eta = aux[0,:] + h
   land = ma.masked_where(h>drytol, eta)
   return land


def depth(current_data):
   """
   Return a masked array containing the depth of fluid only in wet cells.
   """
   from numpy import ma
   drytol = getattr(current_data.user, 'dry_tolerance', drytol_default)
   q = current_data.q
   h = q[0,:]
   depth = ma.masked_where(h<=drytol, h)
   return depth

def surface_full1(current_data):
   #print "full solver"
   """
   Return a masked array containing the surface elevation only in wet cells.
   Surface is eta = h+topo.
   Return the result from the full solver
   """
   from numpy import ma
   drytol = getattr(current_data.user, 'dry_tolerance', drytol_default)
   q = current_data.q
   aux = current_data.aux
   h = q[0,:]
   eta = aux[0,:] + h
   water = ma.masked_where(h<=drytol, eta)
   return water
   
def surface_full2(current_data):
   """
   Return a masked array containing the surface elevation only in wet cells.
   Also masked to show only where the full solver is used
   Surface is eta = h+topo.
   Return the result from the full solver
   """
   from numpy import ma
   drytol=getattr(current_data.user,'dry_tolderance',drytol_default)
   q=current_data.q
   aux=current_data.aux
   h=q[2,:]
   eta=aux[0,:]+h
   water=ma.masked_where(h<=drytol,eta)
   water=ma.masked_where(aux[1,:]!=2.0,water)
   return water
   
def surface_roe2(current_data):
   #print "roe solver"
   """
   Return a masked array containing the surface elevation only in wet cells.
   Also masked to show only where the roe solver is used
   Surface is eta = h+topo.
   Return the result from the full solver
   """
   from numpy import ma
   drytol=getattr(current_data.user,'dry_tolderance',drytol_default)
   q=current_data.q
   aux=current_data.aux
   h=q[2,:]
   eta=aux[0,:]+h
   water=ma.masked_where(h<=drytol,eta)
   water=ma.masked_where(aux[1,:]!=1.0,water)
   return water
   
def roe_error(current_data):
	#print "roe error"
	"""
	Return an array containing the absolute value of the difference
	between the roe result and the full Riemann result
	"""
	import numpy as np
	q=current_data.q
	aux=current_data.aux
	hfull=q[0,:]
	hroe=q[2,:]
	error=abs(hfull-hroe)
	return error
   
   
   
   
   
   
   
   
   
   
   
   


