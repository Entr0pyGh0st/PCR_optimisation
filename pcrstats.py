# Calculates metrics for a finished PCR run
#
# Ben Shirt-Ediss, August 2017
# http://shirt-ediss.me

import pcrparam


def get_state_variable(table, column_num) :
	"""
	If a table is represented as a list of equal sized row lists,
	this extracts a column of that table, as a list
	"""
	return [ row[column_num] for row in table ]




def stats_yield(allstate, allspecies, param) :
	"""
	Yield in ng/ul of newly synthesised dsDNA amplicon
	"""

	amplicon_index = allspecies.index("DNA")

	amplicon_0 = allstate[0][amplicon_index]
	amplicon_final = allstate[-1][amplicon_index]

	amplicon_new_conc_M = amplicon_final - amplicon_0

	if amplicon_new_conc_M > 0 :
		return pcrparam.molar_nM_to_mass_ngul(amplicon_new_conc_M * 1e9, param.amplicon_length, 0)
	else :
		return 0




def stats_times_amplification(allstate, allspecies, param) :

	amplicon_index = allspecies.index("DNA")

	amplicon_0 = allstate[0][amplicon_index]
	amplicon_final = allstate[-1][amplicon_index]

	if amplicon_0 == 0 :
		return 0
	else :
		return amplicon_final / amplicon_0




def stats_purity_dna(allstate, allspecies, param) :
	"""
	Purity of DNA amplicon in final state: mass conc of DNA amplicon / total mass conc of "DNA like" species

	MUST use mass concs here: molar concs don't properly reveal where the material in the system has gone
	"""

	final_state = allstate[-1]

	DNA_index = allspecies.index("DNA")
	DNA_nM_conc = final_state[DNA_index] * 1e9
	DNA_ngul_conc = pcrparam.molar_nM_to_mass_ngul(DNA_nM_conc, param.amplicon_length, 0)

	all_others_ngul_conc = 0
	for i, species in enumerate(allspecies) :
		if species != "dNTP" :	# can add others here, like primers and single strands?
			nt, bp = param.nt_bp_count(species)
			species_nM_conc = final_state[i] * 1e9
			species_ngul_conc = pcrparam.molar_nM_to_mass_ngul(species_nM_conc, bp, nt)

			all_others_ngul_conc += species_ngul_conc

	return DNA_ngul_conc / all_others_ngul_conc




def stats_purity_dnaOLD(allstate, allspecies, param) :
	"""
	Returns absolute purity of amplicon, in terms of percent of nt's that are in amplicons
	But harold said: dNTPs, primers... are not going to be useful in downstream processes
	and so the purity should be the amount of DNA amplicon / total amount of "DNA like" species
	"""

	amplicon_index = allspecies.index("DNA")
	primertop_index = allspecies.index("PT0")
	primerbottom_index = allspecies.index("PB0")
	dntp_index = allspecies.index("dNTP")

	amplicon_0 = allstate[0][amplicon_index]
	amplicon_final = allstate[-1][amplicon_index]

	amplicon_new_conc_M = amplicon_final - amplicon_0

	if amplicon_new_conc_M > 0 :
		# conc of bases in newly synthesised amplicons
		CBamplicon = 2 * param.amplicon_length * amplicon_new_conc_M

		primertop_0 = allstate[0][primertop_index]
		primerbottom_0 = allstate[0][primerbottom_index]
		dntp_0 = allstate[0][dntp_index]

		# total conc of single bases in closed system, if all complexes were broken up into single bases
		# best calculated from initial condition. Initial amplicon is double stranded, primers and dNTP are single stranded
		CBsystem = (2 * param.plasmid_length * amplicon_0) + (param.primer_length * primertop_0) + (param.primer_length * primerbottom_0) + (1 * dntp_0)

		return CBamplicon / CBsystem
	else :
		return 0




def stats_final_state_dict(allstate, allspecies, param) :
	"""
	Makes a dictionary from the final sim state, for gel display
	species: [nt, bp, final molar conc, final bp mass conc]
	"""

	fsd = {}

	for i, species in enumerate(allspecies) :

		nt, bp = param.nt_bp_count(species)
		molar_conc = allstate[-1][i]

		fsd[species] = [nt, bp, molar_conc, pcrparam.molar_nM_to_mass_ngul(molar_conc * 1e9, bp, nt)]

			# note: the final mass conc is just worked out from
			# the dsDNA fraction of a complex -- can be inaccurate, but good enough for DoE workshop

			# Not implemented: the DNA species has initial conc in plasmids taken off
			# so it refers to free-floating DNA

	return fsd

 
























#
#
# EXECUTE
#
#

def analyse_results(alltime, allstate, allspecies, param) :

	trajectory = {}

	# timeseries data: for drawing amplification curve
	trajectory['time'] = alltime
	amplicon_index = allspecies.index("DNA")
	trajectory['target_amplicon_molar'] = get_state_variable(allstate, amplicon_index)

	# final state: for drawing the bands on the gel
	trajectory['final_state'] = stats_final_state_dict(allstate, allspecies, param)

	# amplification stats: for display on interface
	trajectory['statistics'] = {
		'yield' : stats_yield(allstate, allspecies, param),
		'purity' : stats_purity_dna(allstate, allspecies, param),
		'times_amplification' : stats_times_amplification(allstate, allspecies, param),
		'runtime_sec' : int(alltime[-1])
	}

	return trajectory

