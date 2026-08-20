"""




"""

from urllib.request import urlopen, Request
from bs4 import BeautifulSoup
import pandas as pd
import pyDOE
import numpy as np
import csv
import time

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


def dataframe_generate(data, rowlabels: list, columnlabels: list, categorical=False):
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


def dataframe_update_values(factorDF, value_list: list):
    """
    takes value_list and updates the rows of the factorDF DataFrame object using dataframe.itterrows()
    returns an updated factorDF
    """
    for index, row in factorDF.iterrows():
        if index == "value":
            for i in range((len(row))):
                row[i] = value_list[i]
            # for keys in integer_factors:
            # row[integer_factors[keys]] = int(row[integer_factors[keys]])
    return factorDF


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


def sukharev(input_values, nr_factors, base, return_design=True, **kwargs):
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


def sobol(data_package,return_design=True, **kwargs):
    """
    :param *args: [min_values,max_values,runs,nr_factors]
    :param kwargs:
    :return: an updated DoE design and a boilerplate DoE design.
    """

    DOE_design = pyDOE.sobol_sequence(data_package[2], data_package[3])
    DOE_design2 = pyDOE.sobol_sequence(data_package[2], data_package[3])

    for i in range(data_package[3]):
        DOE_design[:, i] *= data_package[1][i]
    if not return_design:
        return DOE_design
    else:
        return DOE_design, DOE_design2


def update_DOEmatrix_datatypes_int64(DOE_matrix, **kwargs):
    """
    takes the DOE_matrix in np.array64, takes a {factor_name:column_index} dictionary, and changes
    the DOE_matrix[:column_index] to data_type ("int","np.int32","np.int64")

        ## NUMPY DOESNT ALLOW DIFFERENT COLUMN FORMATS. HOW DID THIS EVEN WORK BEFORE? ##

    :param DOE_matrix: updated DOE design
    :param data_type:"int" OR "np.int32" OR "np.int64"
    :param kwargs: {factor_name:column_index}
    :return: updated DOE_matrix
    """
    for keys in kwargs:
        for rows in DOE_matrix[:, kwargs[keys]]:
            DOE_matrix[int(rows), kwargs[keys]].astype(np.int64)
    return DOE_matrix

def DOE_results_to_list(DOE_output):
    result_yield = []
    result_purity = []
    result_amplificationpct = []
    result_duration = []

    for b in DOE_output[
        0].keys():  # DOE output is a dictionary, whose keys are the run nr. and the values a dictionary.
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
    data_by_columns = [result_yield, result_purity, result_amplificationpct, result_duration]
    data_by_rows = list(map(list, zip(*data_by_columns)))


def save_data_to_csv(data, filename, factorinfoDF, DOE_matrix):
    csvheader = list(factorinfoDF.columns)
    for row in output_labels:
        csvheader.append(row)

    csvbody = np.concatenate([DOE_matrix, data], axis=1)

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
    x_axis = [str(i) for i in range(len(data[0]))]  ## any result out_file goes really

    axs[0, 0].scatter(x_axis, data[0], marker=".")
    axs[0, 1].scatter(x_axis, data[1])
    axs[1, 0].scatter(x_axis, data[2])
    axs[1, 1].scatter(x_axis, data[3])

    plt.show()


def show_plots_byFactor(DOE_matrix, data, column=1):
    counter = 0

    fg, axs = plt.subplots(nrows=4, ncols=3, figsize=(5.5, 3.5), layout="constrained")

    y_axis = data[column]

    for row in range(len(axs)):
        for column in range(len(axs[row])):
            x_axis = DOE_matrix[:, counter]
            axs[row, column].scatter(x_axis, y_axis)
            axs[row, column].set_title(factorList[counter])
            counter += 1

    plt.show()


# show_plots_outputs()

class DataBall:
    """
    dataBall is a data handler and repository for any DoE campaign.
    it holds:
        - data on starting parameters.
        - data on boilerplate and updated DoE designs.
        - data on simulated DoE runs.
        -
    Instantiates, formats and holds starting parameter information.
    Instantiates, formats and holds DoE designs and parametrized DoE designs (i.e. DoE matrix).

    """

    def __init__(self):
        self.HTMLfile = open("C:/Users/gonca/Desktop/Python/PcrOptimiser/VirtualPCRSimulator.html")
        self.continuousFactors, self.categoricalFactors = BSParse(self.HTMLfile)

        self.factor_names = HTMLDataExtractor(self.continuousFactors, "name")  # extracts names from the website but these are different from the source documentation
        self.factor_count = len(self.factor_names)
        self.factor_attributes = ["min", "max", "value", "type"]
        self.factor_min = HTMLDataExtractor(self.continuousFactors, "min")
        self.factor_max = HTMLDataExtractor(self.continuousFactors, "max")
        self.factor_setvalues = HTMLDataExtractor(self.continuousFactors, "value")
        self.factor_datatype = HTMLDataExtractor(self.continuousFactors, "type")
        self.factor_polymerase = HTMLDataExtractor(self.continuousFactors, "value")
        self.factor_class_integers = {"id_cycles": 0}
        self.datapackage = [self.factor_min, self.factor_max, self.factor_setvalues, self.factor_datatype]

        self.factor_DF = dataframe_generate(self.datapackage, self.factor_attributes, self.factor_names)
        self.factor_names = [factors for factors in self.factor_DF.columns] # updates factor names with the true ones.

        self.DOE_cache = []  # Structure: [ [index:int,"DOE_<version> - <DOE_design>", NxM matrix:list of lists] ]
        self.DOE_version = 0
        self.DOE_active_pointer = 0
    def function_mapping(self, function):
        function_map = \
            {sobol: [self.factor_min, self.factor_max, int(input("How many runs for the sobol?")), self.factor_count],
             sukharev: [],

             }
        return function_map[function]

    def DOE_import(self, design, self_data=True, **kwargs):
        """
        Takes in a function via the design variable.
        Passes the function into function_mapping() to collect the arguments to pass into itself.
        Creates self."design<version>" variable to hold the updated DoE matrix, in int64 format.
        Creates self."design<version>.1" variable to hold the boilerplate DoE matrix.
        Caches the design and primes the DoE into the system (via pointer)
        """

        self.DOE_version += 1
        self.DOE_active_pointer = self.DOE_version

        if self_data:
            data_package = self.function_mapping(design)
            setattr(self, design.__name__ + str(self.DOE_version),  update_DOEmatrix_datatypes_int64(design(data_package)[0], **self.factor_class_integers))
            setattr(self, design.__name__ + str(self.DOE_version) + "a", design(data_package)[1])

        print(self.__class__.__name__,":",design.__name__ + str(self.DOE_version), "generated")
        print(self.__class__.__name__,":",design.__name__ + str(self.DOE_version) + "a","generated ---> blank matrix")

        self.DOE_cache.append([self.DOE_version,design.__name__+str(self.DOE_version), getattr(self, design.__name__ + str(self.DOE_version))])
        self.DOE_active = getattr(self, design.__name__ + str(self.DOE_version))

        if not self_data:
            pass

    def DOE_current_design(self, change=""):
        """
        Returns the DOE designs in cache to the user.
        passing an integer allows the user to choose which DoE design to prime into the system (via pointer).
        :return:
        """
        if change == "":
            for row_nr in range(len(self.DOE_cache)):
                print("DOE nr.: {0} // name: {1}".format(self.DOE_cache[row_nr][0], self.DOE_cache[row_nr][1]))
            print("Active DOE: >>{}<<".format(self.DOE_active_pointer))
            return

        try:
            change_int = int(change)
        except ValueError:
            print("Invalid input: expected an integer or empty string")
            return

        try:
            self.DOE_active = self.DOE_cache[change_int - 1][2] # DOE cache has the matrix in index 2.
            self.DOE_active_pointer = change_int
        except IndexError:
            print("That index is out of bounds")

    def RUN(self, **kwargs):
        """
        Runs Ben's PCR simulator.
        Runs the DOE design assigned by the pointer.

        :param recent_design:
        :param args:
        :return:
        """
        results = []
        try:
            if kwargs["hard_limit"] == int:
                DOE_matrix = self.DOE_active
                DOE_matrix = DOE_matrix[:kwargs["hard_limit"], :]  # truncates the DoE design
        except:
            pass
        for rows in self.DOE_active:
            test_run = dataframe_update_values(self.factor_DF, rows)
            test_run = dataframe_to_pcr_format(test_run)

            results.append(pcrsim.demo(test_run))
            print("{0} out of {1}".format(len(results),len(self.DOE_active)))

        setattr(self,self.DOE_cache[self.DOE_active_pointer][1]+"_data", results)




    def collect(self, function, new_variable_names):
        """
        receives a function and labels the outputs according to pre-defined instructions.
        :param function: any function
        :param kwargs: dict with .values() for labelling the function outputs
        :return: self.values()
        """
        if "dataframe" in function.__name__.lower():
            function.__call__()
            pass


## script execution ----------------------------------

Experimental_design1 = DataBall() # creates the DataBall object

Experimental_design1.DOE_import(sobol) # creates an DOE design
Experimental_design1.RUN() # runs the DOE design



# saves output data into 2 lists: 1 by row (i.e. run nr. in row) 1 by column (i.e. run nr. by column)
output_labels = ["DNA yield (ng/uL)", "DNA Purity (%)", "x amplification", "duration (s)"]
data_by_rows = []
data_by_columns = []
# save_data_to_list(RUN)

# appends output data to bottom of the csv file.
# save_data_to_csv(data_by_rows,"run_results.csv",factorinfoDF,sobol_matrix)

## script execution ----------------------------------

# test_sobol_matrix = pyDOE.sobol_sequence(50,12)
# test_sobol_matrix =sobol(minList,maxList,50,12)
# test_sobol_matrix = update_DOEmatrix_datatypes_int64(test_sobol_matrix,**integer_factors)
# RUN3 = DOE_simulation(test_sobol_matrix,factorinfoDF,hard_limit=3)
# save_data_to_list(RUN3)
# show_plots_byFactor(test_sobol_matrix,data=data_by_columns)
# show_plots_outputs(data=data_by_columns)
