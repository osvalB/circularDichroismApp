import numpy as np
"""

Based on https://www.chem.uci.edu/~dmitryf/manuals/Fundamentals/CD%20practical%20guide.pdf

Absorbance           units are Absorbance
Molar Extinction     units are A*L/mol*cm
Degrees              units are degrees °
Molar ellipticity    units are deg*cm2/dmol
concentration        units are g/l = mg/ml
MolecularWeight      units are g/mol 
Path length of cell  units are cm

"""

# From absorbance to the desired units

def absorbance2milliabsorbance(absorbance,*args):
	return absorbance * 1e3

def absorbance2molarExtinction(absorbance,concentration,pathLength,molecularWeight):
	return absorbance * molecularWeight / (concentration * pathLength)

def absorbance2degrees(absorbance,*args):
	return absorbance * 32.98

def absorbance2millidegrees(absorbance,*args):
	return absorbance * 32980

def absorbance2molarEllipticity(absorbance,concentration,pathLength,molecularWeight):
	return absorbance2molarExtinction(absorbance,concentration,pathLength,molecularWeight) * 3298

# From the desired units to absorbance

def milliabsorbance2absorbance(milliAbsorbance,*args):
	return milliAbsorbance / 1e3

def molarExtinction2absorbance(molarExtinction,concentration,pathLength,molecularWeight):
	return molarExtinction*concentration*pathLength / molecularWeight

def degrees2absorbance(degrees,*args):
	return degrees / 32.98

def millidegrees2absorbance(millidegrees,*args):
	return millidegrees / 32980

def molarEllipticity2absorbance(molarEllipticity,concentration,pathLength,molecularWeight):
	return molarEllipticity * concentration * pathLength / (3298 * molecularWeight)

## General function to convert between different units

def convert2absorbance(value,unitsBegin,concentration=1,pathLength=1,molecularWeight=50,numberOfCroms=1):

	"""

	Based on https://www.chem.uci.edu/~dmitryf/manuals/Fundamentals/CD%20practical%20guide.pdf
	

	"""

	if unitsBegin in ['meanUnitMolarExtinction','meanUnitMolarEllipticity']:

		try: 
			molecularWeight = molecularWeight / numberOfCroms  
			unitsBegin      = unitsBegin.replace("meanUnitM", "m")

		except:

			# return a numpy array full of NaNs with the same shape as value
			return np.full_like(value, np.nan, dtype=np.float64)

	if unitsBegin == "absorbance":
		return value

	if 0 in [concentration,pathLength,molecularWeight] and unitsBegin in ['molarExtinction','molarEllipticity']:
		# return a numpy array full of NaNs with the same shape as value
		return np.full_like(value, np.nan, dtype=np.float64)

	functions = {
	'milliabsorbance'  : milliabsorbance2absorbance,
	'molarExtinction'  : molarExtinction2absorbance,
	'degrees'          : degrees2absorbance,
	'millidegrees'     : millidegrees2absorbance,
	'molarEllipticity' : molarEllipticity2absorbance
	}

	return functions[unitsBegin](value,concentration,pathLength,molecularWeight)

def absorbance2desiredUnits(value,unitsEnd,concentration=1,pathLength=1,molecularWeight=50,numberOfCroms=1):

	if unitsEnd in ['meanUnitMolarExtinction','meanUnitMolarEllipticity']:

		try: 
			molecularWeight = molecularWeight / numberOfCroms  
			unitsEnd        = unitsEnd.replace("meanUnitM", "m")

		except:
			# return a numpy array full of NaNs with the same shape as value
			return np.full_like(value, np.nan, dtype=np.float64)

	if 0 in [concentration,pathLength,molecularWeight] and unitsEnd in ['molarExtinction','molarEllipticity']:
		# return a numpy array full of NaNs with the same shape as value
		return np.full_like(value, np.nan, dtype=np.float64)
 
	if unitsEnd == "absorbance":
		return value

	functions = {
	'milliabsorbance'  : absorbance2milliabsorbance,
	'molarExtinction'  : absorbance2molarExtinction,
	'degrees'          : absorbance2degrees,
	'millidegrees'     : absorbance2millidegrees,
	'molarEllipticity' : absorbance2molarEllipticity
	}

	return functions[unitsEnd](value,concentration,pathLength,molecularWeight)