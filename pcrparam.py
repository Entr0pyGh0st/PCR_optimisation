# PCR simulator parameters
#
# Ben Shirt-Ediss, August 2017
# http://shirt-ediss.me

class PCRparam :

	# ---------------------- PART 1: HARD SET CONSTANTS ------------------------

	########################################
	# Kinetic Rate Constants
	########################################

	# forward bimolecular rate constant for ALL hybridisation
	kh 						= 2.5e5  			# M^-1 s^-1
		# other studies use 1e6 M^-1 s^-1

	# forward bimolecular rate constant for polymerase association to binary complex
	kTaq 					= 1.9 * 1e7  	# M^-1 s^-1

	# quantities for calculating ke (increase curveH to increase ke)
	ke_curveH 				= 1e5			# dimensionless
			# this value is found to work well with the extension phase

	# quantities for calculating kd (increase A to increase kd)
	kd_A 					= 0.0008
	kd_Ea 					= 800		# cal mol^-1

			# these are just set to give decay of approx half of polymerase by end of experiment

			# info here has denaturation activation energies of other proteins
			# at around 100 kcal/mol = 100000 cal mol^-1
			#https://books.google.co.uk/books?id=SkSQNNACcrYC&pg=PA426&lpg=PA426&dq=activation+energy+for+polymerase+denaturation&source=bl&ots=MqzH27MzKZ&sig=a1PSMASygW_5vpOuTjT5pdQxSgo&hl=en&sa=X&ved=0ahUKEwjUjMbl7IHWAhVkLsAKHYHBCToQ6AEIRzAF#v=onepage&q=activation%20energy%20for%20polymerase%20denaturation&f=false
			# !! but this activation energy value means that no polymerase ever decays


	########################################
	# Stock Concentrations
	########################################

	PT0_stock_conc 			= 10e-6		# M
	PB0_stock_conc 			= 10e-6		# M
	dNTP_stock_conc   		= 10e-3		# M

	########################################
	# Volumes
	########################################

	buffer_vol				= 10  # uL
	reaction_volume 		= 50  # uL  -- total reaction volume is always 50uL

	########################################
	# Plasmid, Amplicon and Primer lengths
	########################################
	
	# the amplicon is supplied as part of a plasmid -- 1 amplicon per plasmid assumed
	# the plasmid length is needed, so that we can calculate the initial conc
	# of amplicon, when the initial plasmid weight is known

	# REMEMBER TO RE-RUN g_pcrsim.py AFTER CHANGING ANY OF THESE
	# AS A NEW REACTION MODEL NEEDS BUILDING

	plasmid_length 			= 5000 	# bp
	amplicon_length 		= 1000  # bp
	primer_length 			= 20    # nt

	extension_bp			= 98	# the number of basepairs a primer is extended by in 1 reaction
									# minimum: 1
									# maximum: extension length ( = amplicon_length - primer_length)
									#			-- this is where full length extension happens in just 1 reaction

									# the extension length must be exactly divisible by extension_bp

	########################################
	# Numerical Solver
	########################################

	# list of integrators here: 
	# http://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html
	#
	# Non re-entrant solvers:		vode, zvode, lsoda
	# Re-entrant solvers:			dopri5, dop853
	#
	integrator_name 		= 'lsoda'

	# approx number of steps to *record* per second of sim time
	# (this is not how many steps the solver actually performs per second: the solver determines its step size itself, adaptively)
	steps_record_per_sec 	= 5

	########################################
	# Thermodynamic & Universal Constants
	########################################

	# deltaS for all reversible DNA hybridisation reactions. 
	# Only deltaH changes with DNA duplex length for the reactions.
	deltaS					= -300				# cal mol^-1

	# deltaG for all reversible Taq polymerase association reactions (Phusion assumed to be the same)
	deltaG_polymerase		= -11000			# cal mol^-1  (11 kcal mol^-1)

	# GAS CONSTANT
	R 						= 1.9872036			# cal mol^-1 K^-1

	# AVOGADRO CONSTANT
	NA 						= 6.02214086e23 	# mol^-1

	########################################
	# Single strand concentration used on the NEB Tm calculator
	# to get the melt temperatures returned the NEB_Tm method below
	# (This is needed to calculate deltaH for reactions)
	########################################
	CT_melt_curves	=		500e-9	# 500nM, same as NEB calculator


	# ---------------------- PART 2: USER SET CONSTANTS VIA CONSTRUCTOR ------------------------

	def __init__(self,

		# these parameters are set to be bad, 
		# so that the user has to play around with them in order to make a successful PCR
		# + Too few cycles
		# + Denat temp too low to break complexes
		# + Anneal and extension temperatures the wrong way around
		# + Extension time too short
		# + Too little plasmid for amplification over 30 cycles
		# + Polymerase, primers and DNTPs set very low, so amplification is quickly limited

		thermocycle_Tc      = [80, 68, 54],  	# Tc can be 0 < Tc <= 100
		thermocycle_sec     = [10, 10, 20],    	# s. durations
		thermocycle_repeats = 15,              	# number of times to repeat above cycle
		hold_Tc             = 5,               	# Celsius. temperature of final hold step
		hold_sec            = 600,             	# s. duration of final hold step

		polymerase          = 'Taq',        # 'Taq' or 'Phusion'

		# The initial mass of plasmid in nanograms (no noise added)
		Plasmid_mass_ng		= 50,
		# initial pipetted volumes
		# The user interface supplies these as noisy volumes
		PT0_vol             = 1,		# uL
		PB0_vol             = 1,		# uL
		dNTP_vol            = 0.5*4,	# uL for **all four** dNTPs (1uL per dNTP)
		E_Units				= 1,		# Enzyme units U

	) :

		########################################
		# Validations
		########################################
		assert len(thermocycle_Tc) == len(thermocycle_sec)
		# extension_bp must be exactly dividable into extension length
		assert ((self.amplicon_length - self.primer_length) % self.extension_bp) == 0
		
		self.extension_length = self.amplicon_length - self.primer_length

		########################################
		# Thermocycle
		########################################
		self.thermocycle_Tc      = thermocycle_Tc
		self.thermocycle_sec     = thermocycle_sec
		self.thermocycle_repeats = thermocycle_repeats
		self.hold_Tc             = hold_Tc
		self.hold_sec            = hold_sec

		#########################################
		# Polymerase specific parameters
		#########################################
		self.polymerase = polymerase

		if polymerase == 'Taq' :
			self.ke_curveV = 1
			self.ke_curveM = 3.53
			self.ke_curveS = 0.27

		elif polymerase == 'Phusion' :
			self.ke_curveV = 2.6
			self.ke_curveM = 3.4
			self.ke_curveS = 0.27	

		else :
			raise ValueError("polymerase must be either 'Taq' or 'Phusion'.")

		########################################
		# Initial Concentrations
		########################################
		self.DNA_0 	= mass_to_nanomolar(Plasmid_mass_ng, self.plasmid_length, self.reaction_volume) / 1e9
		self.dNTP_0 = dNTP_vol * self.dNTP_stock_conc / self.reaction_volume
		self.PT0_0  = PT0_vol  * self.PT0_stock_conc  / self.reaction_volume
		self.PB0_0  = PB0_vol  * self.PB0_stock_conc  / self.reaction_volume
		self.E_0    = U_to_nanomolar(E_Units) / 1e9  # A "rough guess" calculation


	#
	# Temperature cycle
	#

	def tcycle(self) :
		"""
		Creates a simple temperature cycle

		Returns 3 lists :
		ts -- start time (in seconds) of temperature phase
		te -- end time (in seconds) of temperature phase
		tc -- temperature (celsius) during temperature phase

		Not a class method, as it depends how the param object is initialised
		"""

		ts = []; te = []; tc = []
		t = 0

		for r in range(0, self.thermocycle_repeats) :
			for i, Tc in enumerate(self.thermocycle_Tc) :
				ts.append(t)
				te.append(t + self.thermocycle_sec[i])
				tc.append(Tc)
				t = t + self.thermocycle_sec[i]

		ts.append(t)
		te.append(t + self.hold_sec)
		tc.append(self.hold_Tc)

		return ts, te, tc


	#
	# Composition of DNA species in model (number of free nt's, number of bp's for each species)
	#

	def nt_bp_count(self, species) :

		"""
		The species in the model are:

		B, T   Bottom and top strands 							ssDNA, nt = amplicon length
		DNA    Amplicon 										dsDNA, bp = amplicon length
		PBn    Primer for bottom strand, extended by n 			ssDNA, nt = primer length + n
		PTn    Primer for top strand, extended by n 			ssDNA, nt = primer length + n
		Xn     Top strand with bottom primer, extended by n 	ss-dsDNA, bp = primer length + n, nt = amplicon length - primer length - n
		XnE    As above, but also with polymerase attached 		ss-dsDNA+protein, bp = primer length + n, nt = amplicon length - primer length - n
		Yn     Bottom strand with bottom primer, extended by n  ss-dsDNA, bp = primer length + n, nt = amplicon length - primer length - n
		YnE    As above, but also with polymerase attached 		ss-dsDNA+protein, bp = primer length + n, nt = amplicon length - primer length - n
		dNTP   A single dNTP 									ssDNA, nt=1
		"""

		nt = 0; bp = 0
		if species in ["B", "T"] :
			nt = self.amplicon_length
			bp = 0
		elif species == "DNA" :
			nt = 0
			bp = self.amplicon_length
		elif ("PB" in species) or ("PT" in species) :
			nt = self.primer_length + int(species[2:])
			bp = 0
		elif ("X" in species) or ("Y" in species) :
			if species[-1] == "E" :
				species = species[:-1]	# knock off the polymerase
			n = int(species[1:])
			nt = self.amplicon_length - self.primer_length - n
			bp = self.primer_length + n
		elif species == "dNTP" :
			nt = 1
			bp = 0

		return nt, bp






#
# Best-fit function for the Tm calculated by NEB calculator, in Taq and Phusion buffers
#

def NEB_Tm(basepairs, polymerase) :
	"""
	Returns Tm melting temperature for a DNA fragment of a number
	of base pairs in The Taq or Phusion NEB polymerase BUFFER

	These curves were made by regression fitting data from the
	NEB Tm calculator (http://tmcalculator.neb.com/#!/), 
	for 50% GC sequences where single strands where present in concentration 500nM.
	"""
	
	param = []

	if polymerase == 'Taq' :
		param = [0.00221, -0.003396, 83.1336]
	elif polymerase == 'Phusion' :
		# SAME AS TAQ FOR NOW (DNA considered to melt the same in Taq and Phusion buffers)
		param = [0.00221, -0.003396, 83.1336]
		# Below params, even though fitting the NEB Tm calc, make DNA melt temperatures very high, and give poor yield for Phusion
		#param = [0.00124, 0.0006915, 102.195]
	
		# Q5 melting -- higher than Taq
		#param = [0.0019577413269071683, -0.00098061260827027362, 91.246229126666975]

	if len(param) == 0 :
		return -1
	else :
		a,b,c = param
		return -1 * (1 / ((a * basepairs) + b)) + c



#
# Concentration conversion methods
#

def mass_to_nanomolar(ng, bp, V_ul) :
	"""
	The nanomolar concentration resulting when ng nanograms of dsDNA (linear or circular)
	of length bp base pairs is put in microlitre volume V_ul
	"""

	# 1 mol of single base pairs weights 660 grams
	mol = (ng / 1e9) / (660 * bp)	# grams / grams per mole of dsDNA molecules = number of moles of molecules
	M = mol / (V_ul / 1e6)
	return M * 1e9




def molar_nM_to_mass_ngul(nM, bp, nt) :
	"""
	Approx nM to ng/ul, for a hybrid ss-dsDNA fragment with nt FREE UNBOUND nucleotides and bp base pairs

	I pull apart the strands of the hybrid complex, and treat it like it was an ssDNA
	Then I use ssDNA nM to ng/ul conversion

	Most of the output of the PCR simulator is hybrid ss-dsDNA complexes, and so this is needed
	"""

	# number of nucleotides, if this complex were unravelled and treated as an ssDNA
	nt_ssDNA = (2*bp) + nt

	frac = ((nt_ssDNA * 308.97) + 18.02) / 1e6

	return nM * frac




def U_to_nanomolar(U) :
	"""
	Guestimate:
	We assume 1.25U per 50uL to correspond to the nM amount of enzyme
	which gives an acceptable amplification curve over 35 cycles,
	which is about 200nM
	"""
	return U * (200 / 1.25)




#
# Gaussian noise on pipetting errors
#

def add_pipetting_errors(vol) :
	"""Add normal distributed noise to vol in (uL)
	
	Noise is normal distributed with a standard deviation depending
	on the volume provided. Gilson pipette specifications state

		P2 - for pipetting 0.2 (+/- 0.024) to 2 uL (+/- 0.03) 
		P10 - for pipetting 1 (+/- 0.025) to 10uL (+/- 0.1) 
		P20 - for pipetting 2 (+/- 0.1) to 20uL (+/- 0.2) 

	For the given volume, the smallest possible pipette is selected
	and the standard deviation is determined by linear interpolation
	between the two specifications.
	"""
	from random import gauss

	if vol>= 0.2 and vol<= 2 :    # P2
		sigma = (0.03-0.024)*(vol-0.2)/(2-0.2)+0.024
	elif vol >= 1 and vol <= 10 : # P10
		sigma = (0.1-0.025)*(vol-1)/(10-1)+0.025
	elif vol >= 2 and vol <= 20 : # P20
		sigma = (0.2-0.1)*(vol-2)/(20-2)+0.1
	else :
		raise ValueError("vol must be between 0.2 and 20")

	return vol + gauss(0,sigma)





