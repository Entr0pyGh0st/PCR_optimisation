"""

context: Fellerman H. et. al created https://virtual-pcr.ico2s.org/pcr/, a website that outputs the yield and the
purity of a virtual Pcr reaction containing 12 modifiable parameters. The study of large parameter sets is
facilitated by experimental designs which use aliasing to confound strictly independent parameters within a single
test run.

goal:

Extract parameter information from https://virtual-pcr.ico2s.org/pcr/

apply DoE methodology.

submit the design for testing to the website.

import the export out_file from the website.

generate optimal predictive model.

forecast factors for testing.

repeat DoE methodology loop


"""

from urllib.request import urlopen, Request
from bs4 import BeautifulSoup
import pandas as pd
import pyDOE
import numpy as np
import csv

from pcrmachine import pcrparam
from pcrmachine import pcrsim
import matplotlib.pyplot as plt

# -----




def HTMLrequest(url: str):
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
    }
    req = Request(url, headers=headers)
    with urlopen(req) as response:
        html_content = response.read().decode("utf-8")
    return html_content


def BSParse(html_file):
    soup = BeautifulSoup(html_file, "html.parser")
    simulatorID = soup.find("div", id="simulator")
    inputTAG = simulatorID.find("input")
    optionTAG = simulatorID.find_all("option")
    inputDATA = inputTAG.find_all("input")

    return [inputDATA, optionTAG]





def HTMLDataExtractor(lst, tag: str):
    """
    Takes a stripped down HTML out_file parsed with beatifulsoup and extracts the values within the HTML tag.

    :param lst: BSParse(html_file)
    :param tag: string. either "min","max","type" or "value"
    :return: simple list with all values associated with a given tag.
    """

    data = []
    integer_factors = ["id_cycles"]
    if tag == "value":
        for i in lst:
            try:
                if i["id"] in integer_factors:
                    data.append(int(i[tag]))
                else:
                    data.append(float(i[tag]))
            except (KeyError, ValueError, TypeError):
                data.append(None)

    elif tag == "min":
        for i in lst:
            try:
                if i["id"] in integer_factors:
                    data.append(int(i[tag]))
                else:
                    data.append(float(i[tag]))
            except (KeyError, ValueError, TypeError):
                data.append(float(0.1))

    elif tag == "max":
        for i in lst:
            try:
                if i["id"] in integer_factors:
                    data.append(int(i[tag]))
                else:
                    data.append(float(i[tag]))
            except (KeyError, ValueError, TypeError):
                data.append(float(180))


    else:
        for i in lst:
            try:
                data.append(i[tag])
            except (KeyError, TypeError):
                data.append(None)

    return data




"""
generates pandas.DataFrame object from data, row labels and column labels OF THEIR WEBSITE, NOT THE GIT FILES.

can be used for continuous factors with min,max,current value input, data type attributes by passing these into rowlabels
can be used for categorical factors by passing factor levels to rowlabels and setting categorical=True.

"""


def generatedataframe(data, rowlabels: list, columnlabels: list, categorical=False):
    if not categorical:
        factorDF = pd.DataFrame(data, index=rowlabels, columns=columnlabels)
        factorDF.drop(factorDF.columns[[-1]], axis=1,
                      inplace=True)  # drops the last column. it's a problem with the HTML parsing.
        factorDF["buffer_vol"] = factorDF.loc[:, "primer_vol"]
        factorDF.rename(columns={"buffer_vol": "topprimer_vol", "primer_vol": "bottomprimer_vol"}, inplace=True)
    elif categorical:
        factorDF = pd.DataFrame(data, columns=columnlabels)
    return factorDF




"""
pcrObject is an object containing all the parameters necessary for running pcrsim.py. Said parameters are below.
Note1: thermocycle_Tc has a list with the denaturation, annealing and extension temperatures respectively.
Note2: thermocycle_sec has a list with the denaturation, annealing and extension times respectively.
Note3: hold_tc and hold_sec are fixed.


        thermocycle_Tc      = [80, 68, 54],  	# Tc can be 0 < Tc <= 100. [0,1,2] = Tdenaturation,Tannealing,Textension
		thermocycle_sec     = [10, 10, 20],    	# s. durations. [0,1,2] = Sdenaturation,Sannealing,Sextension
		 thermocycle_repeats = 15,              	# number of times to repeat above cycle
		X hold_Tc             = 5,               	# Celsius. temperature of final hold step !NOT A FAcTOR
		X hold_sec            = 600,             	# s. duration of final hold step !NOT A FAcTOR
		polymerase          = 'Taq',        # 'Taq' or 'Phusion'

		# The initial mass of plasmid in nanograms (no noise added)
		Plasmid_mass_ng		= 50,
		# initial pipetted volumes
		# The user interface supplies these as noisy volumes
		PT0_vol             = 1,		# uL
		PB0_vol             = 1,		# uL
		dNTP_vol            = 0.5*4,	# uL for **all four** dNTPs (1uL per dNTP)
		E_Units				= 1,		# Enzyme units U
"""


def dataframe_to_pcr_format(df, row=2):
    """
    takes a pd.DataFrame out_file and maps it back into the accepted format of the pcrparam object of pcrparam.py
    format conversion is necessary to seamlessly use pcrsim.py and call the pcr run simulation later on.
    row=2 extracts the base values of factorinfoDF and can be iterated upon to receive the ith row of the DOE design.
    """
    mapping_file = {
        "thermocycle_repeats": df.loc[:, "cycles"][row],
        "thermocycle_Tc": [df.loc[:, "denaturation_temperature"][row], df.loc[:, "annealing_temperature"][row],
                           df.loc[:, "extension_temperature"][row]],
        "thermocycle_sec": [df.loc[:, "denaturation_duration"][row], df.loc[:, "annealing_duration"][row],
                            df.loc[:, "extension_duration"][row]],
        "PT0_vol": df.loc[:, "topprimer_vol"][row],
        "PB0_vol": df.loc[:, "bottomprimer_vol"][row],
        "dNTP_vol": df.loc[:, "dNTP_vol"][row],
        "Plasmid_mass_ng": df.loc[:, "plasmid_mass"][row],
        "E_Units": df.loc[:, "polymerase_vol"][row]}

    pcr_object = pcrparam.PCRparam(**mapping_file)

    return pcr_object


"""    for keys in mapping_file:
        setattr(pcr_object,keys,mapping_file[keys])
"""


# ---------------------------------
# >>>> functionality testing <<<<
#  1. dataframe_to_pcr_format creates a pcrparam.PcRparam class object and feeds in a mapping out_file via **kwargs
#  to instantiate the pcr simulation variables using the DataFrame object generated in generatedataframe() and not the native ones.

# pcrObject = dataframe_to_pcr_format(factorinfoDF)
# B = pcrsim.demo(pcrObject)

# ---------------------------------

def dataframe_update_values(factorDF, value_list: list):
    """
    takes value_list and updates the rows of the factorDF DataFrame object using dataframe.itterrows()
    returns an updated factorDF
    """
    for index, row in factorDF.iterrows():
        if index == "value":
            for i in range((len(row))):
                row[i] = value_list[i]
            for keys in integer_factors:
                row[integer_factors[keys]] = int(row[integer_factors[keys]])
    return factorDF


# -----
# factorinfoDF = dataframe_update_values(factorinfoDF,newvalue_list)
# -----

def dataframe_extract(factorDF, data_attribute: str):
    """
    takes the factorinfoDF dataframe and returns the row labelled data_attribute.
    :param factorDF:  factorinfoDF DataFrame
    :param data_attribute: row label ("min","max,"value","type")
    :return: list of row values
    """
    a = []
    for index, row in factorDF.iterrows():
        if index == data_attribute:
            for i in range((len(row))):
                a.append(row[i])
    return a

def sukharev(input_values, nr_factors, base, return_design=False, **kwargs):
    """
    takes a list of input values, computes a sukharev grid that's n_factors wide and base**n_factors long and returns
    the grid's ith column multiplied my the ith value of input_value.
    return_design=false returns the updated DOE matrix.
    return_design=true returns the updated DOE matrix and the original sukharev grid matrix.
    :param input_values: list of integers
    :param nr_factors: int
    :param base: int
    :param return_design: false or true
    :return: nparray.(base**nr_factors,nr_factors) ONCE or TWICE depending if return_design = false or true
    """
    DOE_design = pyDOE.sukharev_grid(base ** nr_factors, nr_factors)
    DOE_design2 = pyDOE.sukharev_grid(base ** nr_factors, nr_factors)

    for i in range(nr_factors):
        DOE_design[:, i] *= input_values[i]
    if not return_design:
        return DOE_design
    else:
        return DOE_design, DOE_design2



def sobol(min_values,max_values,runs,nr_factors,return_design=False,**kwargs):

    """
    :param input_values:
    :param runs:
    :param nr_factors:
    :param return_design:
    :param kwargs:
    :return:
    """

    DOE_design = pyDOE.sobol_sequence(runs, nr_factors)
    DOE_design2 = pyDOE.sobol_sequence(runs, nr_factors)

    for i in range(nr_factors):
        DOE_design[:, i] *= max_values[i]
    if not return_design:
        return DOE_design
    else:
        return DOE_design, DOE_design2


def update_DOEmatrix_datatypes_int64(DOE_matrix, **kwargs):
    """
    takes the DOE_matrix in np.array64, takes a {factor_name:column_index} dictionary, and changes
    the DOE_matrix[:column_index] to data_type ("int","np.int32","np.int64")

    :param DOE_matrix: updated DOE design
    :param data_type:"int" OR "np.int32" OR "np.int64"
    :param kwargs: {factor_name:column_index}
    :return: updated DOE_matrix
    """
    for keys in kwargs:
        for rows in DOE_matrix[:, kwargs[keys]]:
            DOE_matrix[int(rows), kwargs[keys]].astype(np.int64)
    return DOE_matrix


def DOE_simulation(DOE_matrix, factorDF, hard_limit=None):
    """
    takes the formatted DOE matrix, takes the test values of the ith row of the DOE matrix and updates the factorDF
     dataframe object, converts the data from the dataframe format to pcr format, pushes the data into the pcr simulator.
     if **Kwargs = {"hard_limit": integer}, only integer rows of the DOE matrix will be used.
    :param DOE_matrix: updated DOE_design
    :param factorDF:  factorinfoDF matrix
    :param kwargs:
    :return:
    """
    results = {}
    ticker = 0
    if type(hard_limit) == int:
        DOE_matrix = DOE_matrix[:hard_limit, :]  # truncates the DoE design
    for rows in DOE_matrix:
        test_run = dataframe_update_values(factorDF, rows)
        test_run = dataframe_to_pcr_format(test_run)

        a = pcrsim.demo(test_run)
        results[len(results.keys())] = a
        ticker += 1
        print(ticker, " out of", len((DOE_matrix)))
    return results



def save_data_to_list(DOE_output):
    result_yield = []
    result_purity = []
    result_amplificationpct = []
    result_duration = []

    for b in DOE_output[0].keys(): # DOE output is a dictionary, whose keys are the run nr. and the values a dictionary.
        if b == "yield":
            for run_nr in DOE_output:
                result_yield.append(DOE_output[run_nr][b])
        if b == "purity":
            for run_nr in DOE_output:
                result_purity.append(DOE_output[run_nr][b])
        if b == "times_amplification":
            for run_nr in DOE_output:
                result_amplificationpct.append(DOE_output[run_nr][b])
        if b == "runtime_sec":
            for run_nr in DOE_output:
                result_duration.append(DOE_output[run_nr][b])

    global data_by_rows, data_by_columns
    data_by_columns = [result_yield,result_purity,result_amplificationpct,result_duration]
    data_by_rows = list(map(list,zip(*data_by_columns)))


def save_data_to_csv(data,filename,factorinfoDF,DOE_matrix):

    csvheader = list(factorinfoDF.columns)
    for row in output_labels:
        csvheader.append(row)

    csvbody = np.concatenate([DOE_matrix,data],axis=1)

    out_file = open(filename, "w", newline="")
    writer = csv.writer(out_file)

    in_file = open(filename, "r", newline="")
    reader = csv.reader(in_file)
    old_data = list(reader)
    in_file.close()

    for row in old_data:
        writer.writerow(row)
    writer.writerows(csvheader)
    writer.writerows(csvbody)
    out_file.close()


def show_plots_outputs(data):

    fg, axs = plt.subplots(nrows=2, ncols=2, figsize=(5.5, 3.5), layout="constrained")
    x_axis = [str(i) for i in range(len(data[0]))] ## any result out_file goes really

    axs[0, 0].scatter(x_axis, data[0],marker=".")
    axs[0, 1].scatter(x_axis,data[1])
    axs[1, 0].scatter(x_axis,data[2])
    axs[1, 1].scatter(x_axis, data[3])

    plt.show()

def show_plots_byFactor(DOE_matrix,data,column=1):
    counter = 0

    fg, axs = plt.subplots(nrows=4, ncols=3, figsize=(5.5, 3.5), layout="constrained")

    y_axis = data[column]

    for row in range(len(axs)):
        for column in range(len(axs[row])):
            x_axis = DOE_matrix[:,counter]
            axs[row,column].scatter(x_axis,y_axis)
            axs[row,column].set_title(factorList[counter])
            counter +=1

    plt.show()

#show_plots_outputs()


## script execution ----------------------------------

# Website data extraction: take out factors and the min,max and current values
# url = "https://virtual-pcr.ico2s.org/pcr/"
local = open("C:/Users/gonca/Desktop/Python/PcrOptimiser/VirtualPCRSimulator.html")

# Extracts all factors from the website in HTML:TAG format
# ws = HTMLrequest(url)
continuousFactors, categoricalFactors = BSParse(local)

# Extracts and separates data from each factor into individual Lists
attributeList = ["min", "max", "value", "type"]
factorList = HTMLDataExtractor(continuousFactors, "name")  # extracts names of factors
minList = HTMLDataExtractor(continuousFactors, "min")  # extracts min values of factors
maxList = HTMLDataExtractor(continuousFactors, "max")  # extracts max value of factors
valueList = HTMLDataExtractor(continuousFactors, "value")
typeList = HTMLDataExtractor(continuousFactors, "type")  # extracts value type of factors
polymeraseList = HTMLDataExtractor(categoricalFactors, "value")

# defines which factors are integers.
integer_factors = {"id_cycles": 0}  # maps the position of id_cycles to its index in whatever list above.

# lists all factor data
factorInfo = [minList,maxList,valueList,typeList]

#  Generates a pd.Dataframe for continuous factors and categorical factors
factorinfoDF = generatedataframe(factorInfo, attributeList, factorList)
factorinfoDFcat = generatedataframe(polymeraseList, "", ["id_polymerase"],categorical=True)  # generates working table for cat. factors

# Generates the DoE designs
sukharev_matrix = sukharev(maxList, len(maxList), 2)
sobol_matrix = sobol(minList,maxList,1200,12)

# Updates the data type of the DoE designs.
sukharev_matrix = update_DOEmatrix_datatypes_int64(sukharev_matrix, **integer_factors)
sobol_matrix = update_DOEmatrix_datatypes_int64(sobol_matrix, **integer_factors)

# Runs the DoE simulation.
#RUN = DOE_simulation(sukharev_matrix, factorinfoDF, hard_limit=5)
RUN = DOE_simulation(sobol_matrix, factorinfoDF,hard_limit=3)


# saves output data into 2 lists: 1 by row (i.e. run nr. in row) 1 by column (i.e. run nr. by column)
output_labels = ["DNA yield (ng/uL)","DNA Purity (%)","x amplification","duration (s)"]
data_by_rows = []
data_by_columns = []
save_data_to_list(RUN)

# appends output data to bottom of the csv file.
#save_data_to_csv(data_by_rows,"run_results.csv",factorinfoDF,sobol_matrix)

## script execution ----------------------------------

#test_sobol_matrix = pyDOE.sobol_sequence(50,12)
#test_sobol_matrix =sobol(minList,maxList,50,12)
#test_sobol_matrix = update_DOEmatrix_datatypes_int64(test_sobol_matrix,**integer_factors)
#RUN3 = DOE_simulation(test_sobol_matrix,factorinfoDF,hard_limit=3)
#save_data_to_list(RUN3)
#show_plots_byFactor(test_sobol_matrix,data=data_by_columns)
#show_plots_outputs(data=data_by_columns)




