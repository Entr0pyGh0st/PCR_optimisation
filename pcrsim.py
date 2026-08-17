# Integrate PCR simulator ODE set over temperature schedule
# Primer extension 20 to 1000 bp in increments of 98 bp
#
# Ben Shirt-Ediss, August 2017, Modified December 2024
# http://shirt-ediss.me
import pcrnotebook
import pcrparam
import pcrstats

import numpy
import math
from scipy.integrate import ode










#
#
# ODE SET
#
#

def ode_set(t, state, odeparams) :
	B, DNA, E, ED, PB0, PB196, PB294, PB392, PB490, PB588, PB686, PB784, PB882, PB98, PT0, PT196, PT294, PT392, PT490, PT588, PT686, PT784, PT882, PT98, T, X0, X0E, X196, X196E, X294, X294E, X392, X392E, X490, X490E, X588, X588E, X686, X686E, X784, X784E, X882, X882E, X98, X98E, Y0, Y0E, Y196, Y196E, Y294, Y294E, Y392, Y392E, Y490, Y490E, Y588, Y588E, Y686, Y686E, Y784, Y784E, Y882, Y882E, Y98, Y98E, dNTP = state

	kTaq, kTaqr, kd, ke, kh, khr20, khr118, khr216, khr314, khr412, khr510, khr608, khr706, khr804, khr902, khr1000 = odeparams

	ddt__B =  - (kh * T * B) + (khr1000 * DNA) - (kh * PT0 * B) + (khr20 * Y0) - (kh * PT98 * B) + (khr118 * Y98) - (kh * PT196 * B) + (khr216 * Y196) - (kh * PT294 * B) + (khr314 * Y294) - (kh * PT392 * B) + (khr412 * Y392) - (kh * PT490 * B) + (khr510 * Y490) - (kh * PT588 * B) + (khr608 * Y588) - (kh * PT686 * B) + (khr706 * Y686) - (kh * PT784 * B) + (khr804 * Y784) - (kh * PT882 * B) + (khr902 * Y882)
	ddt__DNA =  + (kh * T * B) - (khr1000 * DNA) + (ke * X882E * (dNTP/4.0)) + (ke * Y882E * (dNTP/4.0))
	ddt__E =  - (kTaq * X0 * E) + (kTaqr * X0E) - (kTaq * Y0 * E) + (kTaqr * Y0E) - (kd * E) - (kTaq * X98 * E) + (kTaqr * X98E) - (kTaq * Y98 * E) + (kTaqr * Y98E) - (kTaq * X196 * E) + (kTaqr * X196E) - (kTaq * Y196 * E) + (kTaqr * Y196E) - (kTaq * X294 * E) + (kTaqr * X294E) - (kTaq * Y294 * E) + (kTaqr * Y294E) - (kTaq * X392 * E) + (kTaqr * X392E) - (kTaq * Y392 * E) + (kTaqr * Y392E) - (kTaq * X490 * E) + (kTaqr * X490E) - (kTaq * Y490 * E) + (kTaqr * Y490E) - (kTaq * X588 * E) + (kTaqr * X588E) - (kTaq * Y588 * E) + (kTaqr * Y588E) - (kTaq * X686 * E) + (kTaqr * X686E) - (kTaq * Y686 * E) + (kTaqr * Y686E) - (kTaq * X784 * E) + (kTaqr * X784E) - (kTaq * Y784 * E) + (kTaqr * Y784E) - (kTaq * X882 * E) + (kTaqr * X882E) - (kTaq * Y882 * E) + (kTaqr * Y882E) + (ke * X882E * (dNTP/4.0)) + (ke * Y882E * (dNTP/4.0))
	ddt__ED =  + (kd * E)
	ddt__PB0 =  - (kh * T * PB0) + (khr20 * X0)
	ddt__PB196 =  - (kh * T * PB196) + (khr216 * X196)
	ddt__PB294 =  - (kh * T * PB294) + (khr314 * X294)
	ddt__PB392 =  - (kh * T * PB392) + (khr412 * X392)
	ddt__PB490 =  - (kh * T * PB490) + (khr510 * X490)
	ddt__PB588 =  - (kh * T * PB588) + (khr608 * X588)
	ddt__PB686 =  - (kh * T * PB686) + (khr706 * X686)
	ddt__PB784 =  - (kh * T * PB784) + (khr804 * X784)
	ddt__PB882 =  - (kh * T * PB882) + (khr902 * X882)
	ddt__PB98 =  - (kh * T * PB98) + (khr118 * X98)
	ddt__PT0 =  - (kh * PT0 * B) + (khr20 * Y0)
	ddt__PT196 =  - (kh * PT196 * B) + (khr216 * Y196)
	ddt__PT294 =  - (kh * PT294 * B) + (khr314 * Y294)
	ddt__PT392 =  - (kh * PT392 * B) + (khr412 * Y392)
	ddt__PT490 =  - (kh * PT490 * B) + (khr510 * Y490)
	ddt__PT588 =  - (kh * PT588 * B) + (khr608 * Y588)
	ddt__PT686 =  - (kh * PT686 * B) + (khr706 * Y686)
	ddt__PT784 =  - (kh * PT784 * B) + (khr804 * Y784)
	ddt__PT882 =  - (kh * PT882 * B) + (khr902 * Y882)
	ddt__PT98 =  - (kh * PT98 * B) + (khr118 * Y98)
	ddt__T =  - (kh * T * B) + (khr1000 * DNA) - (kh * T * PB0) + (khr20 * X0) - (kh * T * PB98) + (khr118 * X98) - (kh * T * PB196) + (khr216 * X196) - (kh * T * PB294) + (khr314 * X294) - (kh * T * PB392) + (khr412 * X392) - (kh * T * PB490) + (khr510 * X490) - (kh * T * PB588) + (khr608 * X588) - (kh * T * PB686) + (khr706 * X686) - (kh * T * PB784) + (khr804 * X784) - (kh * T * PB882) + (khr902 * X882)
	ddt__X0 =  + (kh * T * PB0) - (khr20 * X0) - (kTaq * X0 * E) + (kTaqr * X0E)
	ddt__X0E =  + (kTaq * X0 * E) - (kTaqr * X0E) - (ke * X0E * (dNTP/4.0))
	ddt__X196 =  - (kTaq * X196 * E) + (kTaqr * X196E) + (kh * T * PB196) - (khr216 * X196)
	ddt__X196E =  + (ke * X98E * (dNTP/4.0)) + (kTaq * X196 * E) - (kTaqr * X196E) - (ke * X196E * (dNTP/4.0))
	ddt__X294 =  - (kTaq * X294 * E) + (kTaqr * X294E) + (kh * T * PB294) - (khr314 * X294)
	ddt__X294E =  + (ke * X196E * (dNTP/4.0)) + (kTaq * X294 * E) - (kTaqr * X294E) - (ke * X294E * (dNTP/4.0))
	ddt__X392 =  - (kTaq * X392 * E) + (kTaqr * X392E) + (kh * T * PB392) - (khr412 * X392)
	ddt__X392E =  + (ke * X294E * (dNTP/4.0)) + (kTaq * X392 * E) - (kTaqr * X392E) - (ke * X392E * (dNTP/4.0))
	ddt__X490 =  - (kTaq * X490 * E) + (kTaqr * X490E) + (kh * T * PB490) - (khr510 * X490)
	ddt__X490E =  + (ke * X392E * (dNTP/4.0)) + (kTaq * X490 * E) - (kTaqr * X490E) - (ke * X490E * (dNTP/4.0))
	ddt__X588 =  - (kTaq * X588 * E) + (kTaqr * X588E) + (kh * T * PB588) - (khr608 * X588)
	ddt__X588E =  + (ke * X490E * (dNTP/4.0)) + (kTaq * X588 * E) - (kTaqr * X588E) - (ke * X588E * (dNTP/4.0))
	ddt__X686 =  - (kTaq * X686 * E) + (kTaqr * X686E) + (kh * T * PB686) - (khr706 * X686)
	ddt__X686E =  + (ke * X588E * (dNTP/4.0)) + (kTaq * X686 * E) - (kTaqr * X686E) - (ke * X686E * (dNTP/4.0))
	ddt__X784 =  - (kTaq * X784 * E) + (kTaqr * X784E) + (kh * T * PB784) - (khr804 * X784)
	ddt__X784E =  + (ke * X686E * (dNTP/4.0)) + (kTaq * X784 * E) - (kTaqr * X784E) - (ke * X784E * (dNTP/4.0))
	ddt__X882 =  - (kTaq * X882 * E) + (kTaqr * X882E) + (kh * T * PB882) - (khr902 * X882)
	ddt__X882E =  + (ke * X784E * (dNTP/4.0)) + (kTaq * X882 * E) - (kTaqr * X882E) - (ke * X882E * (dNTP/4.0))
	ddt__X98 =  - (kTaq * X98 * E) + (kTaqr * X98E) + (kh * T * PB98) - (khr118 * X98)
	ddt__X98E =  + (ke * X0E * (dNTP/4.0)) + (kTaq * X98 * E) - (kTaqr * X98E) - (ke * X98E * (dNTP/4.0))
	ddt__Y0 =  + (kh * PT0 * B) - (khr20 * Y0) - (kTaq * Y0 * E) + (kTaqr * Y0E)
	ddt__Y0E =  + (kTaq * Y0 * E) - (kTaqr * Y0E) - (ke * Y0E * (dNTP/4.0))
	ddt__Y196 =  - (kTaq * Y196 * E) + (kTaqr * Y196E) + (kh * PT196 * B) - (khr216 * Y196)
	ddt__Y196E =  + (ke * Y98E * (dNTP/4.0)) + (kTaq * Y196 * E) - (kTaqr * Y196E) - (ke * Y196E * (dNTP/4.0))
	ddt__Y294 =  - (kTaq * Y294 * E) + (kTaqr * Y294E) + (kh * PT294 * B) - (khr314 * Y294)
	ddt__Y294E =  + (ke * Y196E * (dNTP/4.0)) + (kTaq * Y294 * E) - (kTaqr * Y294E) - (ke * Y294E * (dNTP/4.0))
	ddt__Y392 =  - (kTaq * Y392 * E) + (kTaqr * Y392E) + (kh * PT392 * B) - (khr412 * Y392)
	ddt__Y392E =  + (ke * Y294E * (dNTP/4.0)) + (kTaq * Y392 * E) - (kTaqr * Y392E) - (ke * Y392E * (dNTP/4.0))
	ddt__Y490 =  - (kTaq * Y490 * E) + (kTaqr * Y490E) + (kh * PT490 * B) - (khr510 * Y490)
	ddt__Y490E =  + (ke * Y392E * (dNTP/4.0)) + (kTaq * Y490 * E) - (kTaqr * Y490E) - (ke * Y490E * (dNTP/4.0))
	ddt__Y588 =  - (kTaq * Y588 * E) + (kTaqr * Y588E) + (kh * PT588 * B) - (khr608 * Y588)
	ddt__Y588E =  + (ke * Y490E * (dNTP/4.0)) + (kTaq * Y588 * E) - (kTaqr * Y588E) - (ke * Y588E * (dNTP/4.0))
	ddt__Y686 =  - (kTaq * Y686 * E) + (kTaqr * Y686E) + (kh * PT686 * B) - (khr706 * Y686)
	ddt__Y686E =  + (ke * Y588E * (dNTP/4.0)) + (kTaq * Y686 * E) - (kTaqr * Y686E) - (ke * Y686E * (dNTP/4.0))
	ddt__Y784 =  - (kTaq * Y784 * E) + (kTaqr * Y784E) + (kh * PT784 * B) - (khr804 * Y784)
	ddt__Y784E =  + (ke * Y686E * (dNTP/4.0)) + (kTaq * Y784 * E) - (kTaqr * Y784E) - (ke * Y784E * (dNTP/4.0))
	ddt__Y882 =  - (kTaq * Y882 * E) + (kTaqr * Y882E) + (kh * PT882 * B) - (khr902 * Y882)
	ddt__Y882E =  + (ke * Y784E * (dNTP/4.0)) + (kTaq * Y882 * E) - (kTaqr * Y882E) - (ke * Y882E * (dNTP/4.0))
	ddt__Y98 =  - (kTaq * Y98 * E) + (kTaqr * Y98E) + (kh * PT98 * B) - (khr118 * Y98)
	ddt__Y98E =  + (ke * Y0E * (dNTP/4.0)) + (kTaq * Y98 * E) - (kTaqr * Y98E) - (ke * Y98E * (dNTP/4.0))
	ddt__dNTP =  - (98 * (ke * X0E * (dNTP/4.0))) - (98 * (ke * Y0E * (dNTP/4.0))) - (98 * (ke * X98E * (dNTP/4.0))) - (98 * (ke * Y98E * (dNTP/4.0))) - (98 * (ke * X196E * (dNTP/4.0))) - (98 * (ke * Y196E * (dNTP/4.0))) - (98 * (ke * X294E * (dNTP/4.0))) - (98 * (ke * Y294E * (dNTP/4.0))) - (98 * (ke * X392E * (dNTP/4.0))) - (98 * (ke * Y392E * (dNTP/4.0))) - (98 * (ke * X490E * (dNTP/4.0))) - (98 * (ke * Y490E * (dNTP/4.0))) - (98 * (ke * X588E * (dNTP/4.0))) - (98 * (ke * Y588E * (dNTP/4.0))) - (98 * (ke * X686E * (dNTP/4.0))) - (98 * (ke * Y686E * (dNTP/4.0))) - (98 * (ke * X784E * (dNTP/4.0))) - (98 * (ke * Y784E * (dNTP/4.0))) - (98 * (ke * X882E * (dNTP/4.0))) - (98 * (ke * Y882E * (dNTP/4.0)))


	return [ddt__B, ddt__DNA, ddt__E, ddt__ED, ddt__PB0, ddt__PB196, ddt__PB294, ddt__PB392, ddt__PB490, ddt__PB588, ddt__PB686, ddt__PB784, ddt__PB882, ddt__PB98, ddt__PT0, ddt__PT196, ddt__PT294, ddt__PT392, ddt__PT490, ddt__PT588, ddt__PT686, ddt__PT784, ddt__PT882, ddt__PT98, ddt__T, ddt__X0, ddt__X0E, ddt__X196, ddt__X196E, ddt__X294, ddt__X294E, ddt__X392, ddt__X392E, ddt__X490, ddt__X490E, ddt__X588, ddt__X588E, ddt__X686, ddt__X686E, ddt__X784, ddt__X784E, ddt__X882, ddt__X882E, ddt__X98, ddt__X98E, ddt__Y0, ddt__Y0E, ddt__Y196, ddt__Y196E, ddt__Y294, ddt__Y294E, ddt__Y392, ddt__Y392E, ddt__Y490, ddt__Y490E, ddt__Y588, ddt__Y588E, ddt__Y686, ddt__Y686E, ddt__Y784, ddt__Y784E, ddt__Y882, ddt__Y882E, ddt__Y98, ddt__Y98E, ddt__dNTP]




def integrate(ode_set, state0, time, odeparams, param) :

	yy = []
	yy.append(state0)
	t0 = time[0]

	solver = ode(ode_set)
	solver.set_integrator(param.integrator_name)
	solver.set_initial_value(state0, t0)
	solver.set_f_params(odeparams)
	for t in time[1:] :
		yy.append(solver.integrate(t))
		if not solver.successful() :
			print("Warning. There are some integration errors.")

	return yy









































#
#
# RATE CONSTANT CALCULATION
#
#

def Keq_vant_hoff(Tc, deltaH, deltaS, R) :

	# Returns equilibrium constant for a reversible reaction, 
	# given the thermodynamic parameters of the reaction are deltaH and deltaS

	# Tc: 			Celsius
	# deltaH:		cal mol^-1
	# deltaS:		cal mol^-1
	# R:			cal mol^-1 K^-1
	
	T = Tc + 273.15

	frac1 = -1 * (deltaH / (R * T))
	frac2 = deltaS / R

	return math.exp(frac1 + frac2)




def calc_hybridisation_reverse_rate_constant(kf, Tc, deltaH, deltaS, R) :

	# Returns the reverse rate constant for a hybridisation reaction (unit: s^-1),
	# given the forward rate constant (and Keq parameters)

	# kf:			M^-1 s^-1
	# Tc:			Celsius
	# deltaH:		cal mol^-1
	# deltaS:		cal mol^-1
	# R:			cal mol^-1 K^-1

	return kf / Keq_vant_hoff(Tc, deltaH, deltaS, R)




def calc_deltaH(Tm, deltaS, R, CT) :
	
	# Given the melting temperature (in celsius) for a hybridisation reaction
	# (measured at total single strand molar concentration CT),
	# and the deltaS for that reaction, returns what the deltaH must be
	#
	# NOTE: deltaH is returned in cal mol^-1

	# Tm:			Celsius
	# deltaS:		cal mol^-1
	# R:			cal mol^-1 K^-1
	# CT:			M
	
	T = Tm + 273.15

	return T * (deltaS - (R * math.log(4.0/CT)))




def calc_polymerase_reverse_rate_constant(kf, Tc, deltaG, R) :
	
	# Calc polymerase off-rate, given polymer on-rate and binding thermodynamics
	
	# kf:			M^-1 s^-1
	# Tc:			Celsius
	# deltaG:		cal mol^-1
	# R: 			cal mol^-1 K^-1

	T = Tc + 273.15
	
	Keq = math.exp((-1 * deltaG) / (R * T))

	return kf / Keq




def lognormal(Tc, M, S) :
	
	# The log-normal curve. The horizontally flipped version of this 
	# is used to work out polymerase extension kinetic forward rate, given temperature
	
	# Tc:		Celsius

	epower = (-1 * (math.log(Tc) - M)**2) / (2 * S**2)
	coeff = 1 / (S * math.sqrt(2 * math.pi) * Tc)

	return coeff * math.exp(epower)




def calc_ke(Tc, param) :

	# Tc:		Celsius
	
	return param.ke_curveH * param.ke_curveV * lognormal((100.0001-Tc), param.ke_curveM, param.ke_curveS)
	



def calc_kd(Tc, param) :

	# Tc:		Celsius

	T = Tc + 273.15

	return param.kd_A * math.exp((-1 * param.kd_Ea) / (param.R * T))









































#
#
# RUN SIMULATION
#
#

def runpcr(param) :

	alltime = []
	allstate = []
	allTc = []
	allcycleinfo = []

	# (1) INITIAL STATE (Molar)

	state0 = [
		0, 				 # B
		param.DNA_0, 			 # DNA
		param.E_0, 			 # E
		0, 				 # ED
		param.PB0_0, 			 # PB0
		0, 				 # PB196
		0, 				 # PB294
		0, 				 # PB392
		0, 				 # PB490
		0, 				 # PB588
		0, 				 # PB686
		0, 				 # PB784
		0, 				 # PB882
		0, 				 # PB98
		param.PT0_0, 			 # PT0
		0, 				 # PT196
		0, 				 # PT294
		0, 				 # PT392
		0, 				 # PT490
		0, 				 # PT588
		0, 				 # PT686
		0, 				 # PT784
		0, 				 # PT882
		0, 				 # PT98
		0, 				 # T
		0, 				 # X0
		0, 				 # X0E
		0, 				 # X196
		0, 				 # X196E
		0, 				 # X294
		0, 				 # X294E
		0, 				 # X392
		0, 				 # X392E
		0, 				 # X490
		0, 				 # X490E
		0, 				 # X588
		0, 				 # X588E
		0, 				 # X686
		0, 				 # X686E
		0, 				 # X784
		0, 				 # X784E
		0, 				 # X882
		0, 				 # X882E
		0, 				 # X98
		0, 				 # X98E
		0, 				 # Y0
		0, 				 # Y0E
		0, 				 # Y196
		0, 				 # Y196E
		0, 				 # Y294
		0, 				 # Y294E
		0, 				 # Y392
		0, 				 # Y392E
		0, 				 # Y490
		0, 				 # Y490E
		0, 				 # Y588
		0, 				 # Y588E
		0, 				 # Y686
		0, 				 # Y686E
		0, 				 # Y784
		0, 				 # Y784E
		0, 				 # Y882
		0, 				 # Y882E
		0, 				 # Y98
		0, 				 # Y98E
		param.dNTP_0] 			 # dNTP

	allspecies = ["B","DNA","E","ED","PB0","PB196","PB294","PB392","PB490","PB588","PB686","PB784","PB882","PB98","PT0","PT196","PT294","PT392","PT490","PT588","PT686","PT784","PT882","PT98","T","X0","X0E","X196","X196E","X294","X294E","X392","X392E","X490","X490E","X588","X588E","X686","X686E","X784","X784E","X882","X882E","X98","X98E","Y0","Y0E","Y196","Y196E","Y294","Y294E","Y392","Y392E","Y490","Y490E","Y588","Y588E","Y686","Y686E","Y784","Y784E","Y882","Y882E","Y98","Y98E","dNTP"]
	allcycleinfo_cols = ["start", "end", "Tc", "kTaq", "kTaqr", "kd", "ke", "kh", "khr20", "khr118", "khr216", "khr314", "khr412", "khr510", "khr608", "khr706", "khr804", "khr902", "khr1000"]
	
	# (2) MAKE TEMPERATURE SCHEDULE LIST

	ts, te, tc = param.tcycle()

	# (3) GO! 

	state = state0

	for i, Tc in enumerate(tc) :

		# CALC RATE CONSTANTS FOR CURRENT TEMPERATURE
		kTaq = param.kTaq
		kTaqr = calc_polymerase_reverse_rate_constant(param.kTaq, Tc, param.deltaG_polymerase, param.R)
		kd = calc_kd(Tc, param)
		ke = calc_ke(Tc, param)
		kh = param.kh

		Tm = pcrparam.NEB_Tm(20, param.polymerase)
		deltaH = calc_deltaH(Tm, param.deltaS, param.R, param.CT_melt_curves)
		khr20 = calc_hybridisation_reverse_rate_constant(param.kh, Tc, deltaH, param.deltaS, param.R)
		Tm = pcrparam.NEB_Tm(118, param.polymerase)
		deltaH = calc_deltaH(Tm, param.deltaS, param.R, param.CT_melt_curves)
		khr118 = calc_hybridisation_reverse_rate_constant(param.kh, Tc, deltaH, param.deltaS, param.R)
		Tm = pcrparam.NEB_Tm(216, param.polymerase)
		deltaH = calc_deltaH(Tm, param.deltaS, param.R, param.CT_melt_curves)
		khr216 = calc_hybridisation_reverse_rate_constant(param.kh, Tc, deltaH, param.deltaS, param.R)
		Tm = pcrparam.NEB_Tm(314, param.polymerase)
		deltaH = calc_deltaH(Tm, param.deltaS, param.R, param.CT_melt_curves)
		khr314 = calc_hybridisation_reverse_rate_constant(param.kh, Tc, deltaH, param.deltaS, param.R)
		Tm = pcrparam.NEB_Tm(412, param.polymerase)
		deltaH = calc_deltaH(Tm, param.deltaS, param.R, param.CT_melt_curves)
		khr412 = calc_hybridisation_reverse_rate_constant(param.kh, Tc, deltaH, param.deltaS, param.R)
		Tm = pcrparam.NEB_Tm(510, param.polymerase)
		deltaH = calc_deltaH(Tm, param.deltaS, param.R, param.CT_melt_curves)
		khr510 = calc_hybridisation_reverse_rate_constant(param.kh, Tc, deltaH, param.deltaS, param.R)
		Tm = pcrparam.NEB_Tm(608, param.polymerase)
		deltaH = calc_deltaH(Tm, param.deltaS, param.R, param.CT_melt_curves)
		khr608 = calc_hybridisation_reverse_rate_constant(param.kh, Tc, deltaH, param.deltaS, param.R)
		Tm = pcrparam.NEB_Tm(706, param.polymerase)
		deltaH = calc_deltaH(Tm, param.deltaS, param.R, param.CT_melt_curves)
		khr706 = calc_hybridisation_reverse_rate_constant(param.kh, Tc, deltaH, param.deltaS, param.R)
		Tm = pcrparam.NEB_Tm(804, param.polymerase)
		deltaH = calc_deltaH(Tm, param.deltaS, param.R, param.CT_melt_curves)
		khr804 = calc_hybridisation_reverse_rate_constant(param.kh, Tc, deltaH, param.deltaS, param.R)
		Tm = pcrparam.NEB_Tm(902, param.polymerase)
		deltaH = calc_deltaH(Tm, param.deltaS, param.R, param.CT_melt_curves)
		khr902 = calc_hybridisation_reverse_rate_constant(param.kh, Tc, deltaH, param.deltaS, param.R)
		Tm = pcrparam.NEB_Tm(1000, param.polymerase)
		deltaH = calc_deltaH(Tm, param.deltaS, param.R, param.CT_melt_curves)
		khr1000 = calc_hybridisation_reverse_rate_constant(param.kh, Tc, deltaH, param.deltaS, param.R)

		odeparams = [kTaq, kTaqr, kd, ke, kh, khr20, khr118, khr216, khr314, khr412, khr510, khr608, khr706, khr804, khr902, khr1000]

		# TIME RANGE
		interval = te[i] - ts[i]
		num_steps = int(interval * param.steps_record_per_sec)
		if num_steps < 3 : num_steps = 3
		time = numpy.linspace(ts[i], te[i], num_steps)

		# INTEGRATE small temperature segment
		yy = integrate(ode_set, state, time, odeparams, param)

		alltime = alltime + list(time)
		allstate = allstate + list(yy)
		allTc = allTc + (len(time) * [Tc])
		allcycleinfo = allcycleinfo + [[ts[i], te[i], Tc] + odeparams]

		# INITIAL STATE FOR NEXT TEMPERATURE PHASE IS FINAL STATE OF THIS TEMPERATURE PHASE
		state = yy[-1]

	return pcrstats.analyse_results(alltime, allstate, allspecies, param)
 










































#
#
# DEMO
#
#

def demo() :

	# use default parameters in constructor of PCRparam object
	param = pcrparam.PCRparam()

	print('Initial Concs ---------------------------------------')
	print('DNA (nM)   :', param.DNA_0 * 1e9)
	print('PT0 (nM)    :', param.PT0_0 * 1e9)
	print('PB0 (nM)    :', param.PB0_0 * 1e9)
	print('dNTP (nM)  :', param.dNTP_0 * 1e9)
	print('E (nM)     :', param.E_0    * 1e9)
	
	print('')
	print('Running simulation...')

	traj = runpcr(param)

	print('Done.')
	print('')

	print('PCR Results ---------------------------------------')
	print('Yield (ng/ul):                     ', traj['statistics']['yield'])
	print('Purity (%):                        ', traj['statistics']['purity'])
	print('x Amplification:                   ', traj['statistics']['times_amplification'])
	print('Experiment duration (s):           ', traj['statistics']['runtime_sec'])
	print('---------------------------------------')


if __name__ == '__main__' :
	demo()
