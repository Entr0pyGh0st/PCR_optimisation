# Makes PCR interface in Jupyter Notebook and calls PCR simulation scripts
#
# Ben Shirt-Ediss, 2024
# http://shirt-ediss.me

import ipywidgets as widgets
from IPython.display import Image, display, HTML

import pcrparam
import pcrsim

















#
#
# CALCULATE MASTER MIX VOLUME
#
#

def calc_mastermix_components_volume(dntp_slider, primer_slider, plasmid_slider, polymerase_slider, vol_warning_label) :
	# note: the mastermix is always 50uL, it is topped up with water once the components are all added
	# note: this means that the components themselves cannot exceed 50uL in volume

	dNTP_vol       	=   4 		* dntp_slider.value
	primer_vol     	=   2 		* primer_slider.value
	plasmid_vol    	= 	0.03 	* plasmid_slider.value
	polymerase_vol 	=  	0.2 	* polymerase_slider.value

	total_vol 		=	dNTP_vol + primer_vol + plasmid_vol + polymerase_vol

	if total_vol > 50 :
		total_vol = 50
		vol_warning_label.value = "<span style='color: red;'>&#9888; 50uL exceeded: PCR simulation cannot be started.</span>"
	else :
		vol_warning_label.value = ""

	return total_vol































#
#
# RUN PCR BUTTON AND ON CLICK EVENT HANDLER
#
#

def create_run_pcr_button(all_controls) :

	cycles_slider, denatC_slider, denatS_slider, annealC_slider, annealS_slider, extendC_slider, extendS_slider, \
	dntp_slider, primer_slider, plasmid_slider, polymerase_slider, \
	vol_warning_label, yield_text, purity_text, ttime_text = all_controls

	run_button = widgets.Button(
		description		= "Run PCR",
		button_style	= "danger",
		layout 			= widgets.Layout(margin="20px 0px", width="300px")
	)

	status = widgets.Output()

	# on-click handler defined inside this function. Runs pcr sim and reports status message
	def on_button_click(b):

		with status:

			if vol_warning_label.value != "" :
				return  	# there is an excessive master mix volume warning

			# PCR SIMULATION =====================================================
			# 1) clear current results
			yield_text.value = "--"
			purity_text.value = "--"
			ttime_text.value = "--"

			# 2) disable the Run PCR button, so it cannot be pressed again
			run_button.disabled = True

			# 3) make class of PCR parameters (adding pipetting errors to dNTPs and primers)
			pcr_params = pcrparam.PCRparam(
					thermocycle_Tc = [
						denatC_slider.value,
						annealC_slider.value,
						extendC_slider.value,
					],
					thermocycle_sec = [
						denatS_slider.value,
						annealS_slider.value,
						extendS_slider.value,
					],
					thermocycle_repeats = cycles_slider.value,

					dNTP_vol        = 4 * pcrparam.add_pipetting_errors(dntp_slider.value),
					PT0_vol         = pcrparam.add_pipetting_errors(primer_slider.value),
					PB0_vol         = pcrparam.add_pipetting_errors(primer_slider.value),
					Plasmid_mass_ng = plasmid_slider.value,
					E_Units         = polymerase_slider.value,

					polymerase = "Taq",
				)

			# 4) run PCR
			print("**************************************")
			print("  Running simulation. Please wait...  ")
			print("**************************************")
			print(" ")

			result = pcrsim.runpcr(pcr_params)
			
			# 5) display results
			yield_text.value 			= "%.2f" % result['statistics']['yield']
			purity_text.value 			= "%.2f" % (result['statistics']['purity'] * 100)
			ttime_text.value 			= "%.1f" % (result['statistics']['runtime_sec']	/ 60)
			status.clear_output()

			# 6) ready for next PCR simulation
			run_button.disabled = False
			# ====================================================================

	run_button.on_click(on_button_click)

	return run_button, status






























#
#
# PCR CONTROL WIDGETS
#
#

def display_pcr_interface() :

	display(Image(filename='pcr.png', width=250))

	# ------------------------------
	# Thermocycler
	# ------------------------------

	# CONTROLS

	cycles_slider = widgets.IntSlider(
		value		= 15,
		min			= 1,
		max 		= 100,
		step 		= 1,
		description = "Cycles:",
	)

	denatC_slider = widgets.IntSlider(
		value		= 80,
		min			= 1,
		max	 		= 100,
		step 		= 1,
		description = "Denaturing:",
	)

	denatS_slider = widgets.IntSlider(
		value 		= 10,
		min 		= 1,
		max 		= 120,
		step	  	= 1,
		description = "for",
	)

	annealC_slider = widgets.IntSlider(
		value		= 68,
		min			= 1,
		max 		= 100,
		step 		= 1,
		description = "Annealing:",
	)

	annealS_slider = widgets.IntSlider(
		value 		= 10,
		min 		= 1,
		max 		= 120,
		step 		= 1,
		description	= "for",
	)

	extendC_slider = widgets.IntSlider(
		value		= 54,
		min			= 1,
		max			= 100,
		step		= 1,
		description = "Extension:",
	)

	extendS_slider = widgets.IntSlider(
		value		= 20,
		min			= 1,
		max 		= 120,
		step 		= 1,
		description = "for",
	)

	# LABELS

	degc_label = widgets.Label(value="°C")
	sec_label = widgets.Label(value="seconds")

	# DISPLAY

	display( widgets.HTML(value="<b>Thermocycler Settings</b>") )

	display(cycles_slider)

	denatC = widgets.HBox([denatC_slider, degc_label])
	denatS = widgets.HBox([denatS_slider, sec_label])
	row1 = widgets.HBox([denatC, denatS])
	display(row1)

	annealC = widgets.HBox([annealC_slider, degc_label])
	annealS = widgets.HBox([annealS_slider, sec_label])
	row2 = widgets.HBox([annealC, annealS])
	display(row2)

	extendC = widgets.HBox([extendC_slider, degc_label])
	extendS = widgets.HBox([extendS_slider, sec_label])
	row3 = widgets.HBox([extendC, extendS])
	display(row3)

	# ------------------------------
	# Reaction Master Mix
	# ------------------------------

	# CONTROLS

	dntp_slider = widgets.FloatSlider(
		value 		= 0.5,
		min			= 0.2,
		max			= 20.0,
		step		= 0.1,
		description = "Each DNTP:",
		readout  	= True,
		readout_format = ".1f"
	)

	primer_slider = widgets.FloatSlider(
		value		= 1.0,
		min			= 0.2,
		max			= 20.0,
		step		= 0.1,
		description = "Each Primer:",
		readout  	= True,
		readout_format = ".1f"
	)

	plasmid_slider = widgets.FloatSlider(
		value		= 50,
		min	 		= 0,
		max 		= 100,
		step		= 1,
		description = "Plasmid:",
		readout  	= True,
		readout_format = ".1f"
	)

	polymerase_slider = widgets.FloatSlider(
		value 		= 1,
		min 		= 0.2,
		max 		= 20.0,
		step 		= 0.1,
		description = "Polymerase:",
		readout  	= True,
		readout_format = ".1f"
	)

	# LABELS

	microlitre1_label = widgets.Label(value="uL (from 10mM stock)")
	microlitre2_label = widgets.Label(value="uL (from 10uM stock)")
	U_label = widgets.Label(value="U of Taq Polymerase")
	ng_label = widgets.Label(value="ng")

	# DISPLAY

	dntp = widgets.HBox([dntp_slider, microlitre1_label])
	primer = widgets.HBox([primer_slider, microlitre2_label])
	plasmid = widgets.HBox([plasmid_slider, ng_label])
	polymerase = widgets.HBox([polymerase_slider, U_label])

	mix_col1_title = widgets.HTML(value="<b>Reaction Master Mix Components</b>")

	mix_col1 = widgets.VBox([mix_col1_title, dntp, primer, plasmid, polymerase])

	# ------------------------------
	# Total Reaction Volume Check
	# ------------------------------

	# CONTROLS

	vol_slider = widgets.FloatSlider(
					value		= 0.0,
					min			= 0.0,
					max			= 50.0,
					step		= 0.1,
					description = "",
					disabled	= True,  # Disable the slider (user cannot interact)
					layout		= widgets.Layout(width='300px')
	)

	vol_warning_label = widgets.HTML(
					value		= ""
	)

	# set reading on total volume slider...
	vol_slider.value = calc_mastermix_components_volume(dntp_slider, primer_slider, plasmid_slider, polymerase_slider, vol_warning_label)

	# .. and update this reading when component levels in the master mix change
	def update_vol_slider(change):
		vol_slider.value = calc_mastermix_components_volume(dntp_slider, primer_slider, plasmid_slider, polymerase_slider, vol_warning_label)

	dntp_slider.observe(update_vol_slider)
	primer_slider.observe(update_vol_slider)
	plasmid_slider.observe(update_vol_slider)
	polymerase_slider.observe(update_vol_slider)

	# LABELS

	microlitre3_label = widgets.Label(value="uL")

	# DISPLAY

	mix_col2_title = widgets.HTML(value="<b>Total Volume of Master Mix Components (Max 50ul)</b>")

	vol = widgets.HBox([vol_slider, microlitre3_label])

	mix_col2 = widgets.VBox(
					[mix_col2_title, vol, vol_warning_label],
					layout=widgets.Layout(
						border='solid 1px black',
						padding='10px',
						margin='20px'))

	mix = widgets.HBox([mix_col1, mix_col2])

	display(mix)

	# ------------------------------
	# Run PCR Button and Results
	# ------------------------------

	# CONTROLS

	yield_text = widgets.Text(
		value		= "--",
		disabled	= True,
		layout		= widgets.Layout(width='100px', font_size='20px')
	)

	purity_text = widgets.Text(
		value		= "--",
		disabled	= True,
		layout		= widgets.Layout(width='100px', font_size='20px')
	)

	ttime_text = widgets.Text(
		value		= "--",
		disabled	= True,
		layout		= widgets.Layout(width='100px', font_size='20px')
	)

	all_controls = [cycles_slider, denatC_slider, denatS_slider, annealC_slider, annealS_slider, extendC_slider, extendS_slider,
					dntp_slider, primer_slider, plasmid_slider, polymerase_slider,
					vol_warning_label, yield_text, purity_text, ttime_text]

	run_button, status = create_run_pcr_button(all_controls)	

	# LABELS

	yield_label = widgets.Label(value="Yield (ng/uL)")
	purity_label = widgets.Label(value="Purity (%)")
	ttime_label = widgets.Label(value="PCR Duration (minutes)")

	# DISPLAY

	display(run_button, status)

	display( widgets.HTML(value="<b>Results</b>") )

	yieldd = widgets.HBox([yield_text, yield_label])
	display(yieldd)
	
	purity = widgets.HBox([purity_text, purity_label])
	display(purity)

	ttime = widgets.HBox([ttime_text, ttime_label])
	display(ttime)







