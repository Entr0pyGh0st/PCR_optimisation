import string
from urllib.request import urlopen, Request

import dmba
from bs4 import BeautifulSoup
import pandas as pd
import pyDOE
import numpy as np
import csv
import time
import os
import itertools
import scipy.stats as st

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


def dataframe_to_pcr_format(row):
    """
    takes a pd.DataFrame out_file and maps it back into the accepted format of the pcrparam object of pcrparam.py
    format conversion is necessary to seamlessly use pcrsim.py and call the pcr run simulation later on.
    row=2 extracts the base values of factorinfoDF and can be iterated upon to receive the ith row of the DOE design.
    """
    mapping_file = {
        "thermocycle_repeats": row[0],
        "thermocycle_Tc": [row[1], row[3], row[5]],
        "thermocycle_sec": [row[2], row[4], row[6]],
        "PT0_vol": row[7],
        "PB0_vol": row[9],
        "dNTP_vol": row[8],
        "Plasmid_mass_ng": row[10],
        "E_Units": row[11]}

    pcr_object = pcrparam.PCRparam(**mapping_file)

    return pcr_object


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


def sobol(data_package, return_design=True, **kwargs):
    """
    :param *args: [min_values,max_values,runs,nr_factors]
    :param kwargs:
    :return: an updated DoE design and a boilerplate DoE design.

    [self.factor_DF.loc["min"], self.factor_DF.loc["max"], int(input("How many runs for the sobol?")), self.factor_DF.shape[1]]
    """
    DOE_design = pyDOE.sobol_sequence(data_package[2], data_package[3])
    DOE_design2 = pyDOE.sobol_sequence(data_package[2], data_package[3])

    for factor_nr in range(data_package[3]):
        for row_nr in range(len(DOE_design)):
            DOE_design[row_nr, factor_nr] = data_package[0].iloc[factor_nr] + (
                    data_package[1].iloc[factor_nr] - data_package[0].iloc[factor_nr]) * DOE_design[
                                                row_nr, factor_nr]
    if not return_design:
        return DOE_design
    else:
        return DOE_design, DOE_design2


def ff2n(data_package, return_design=True, **kwargs):
    DOE_design = pyDOE.ff2n(data_package[0])
    DOE_design2 = pyDOE.ff2n(data_package[0])

    for i in range(data_package[0]):
        DOE_design[:, i] = DOE_design[:, i] * np.std([data_package[1][i], data_package[2][i]]) + np.average(
            [data_package[1][i], data_package[2][i]])

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
    print("Imported file shape: ", DOE_matrix.shape)

    for keys in kwargs:
        for rows in DOE_matrix:
            rows[kwargs[keys]].astype(np.int64)
    return DOE_matrix


def fractionalFactorial(data_package):  # factors, resolution, min, max
    design1 = pyDOE.fracfact_by_res(data_package[0], data_package[1])
    design2 = pyDOE.fracfact_by_res(data_package[0], data_package[1])

    for factor_nr in range(data_package[0]):
        for row_nr in range(len(design1)):
            design1[row_nr, factor_nr] = data_package[2].iloc[factor_nr] + (
                    data_package[3].iloc[factor_nr] - data_package[2].iloc[factor_nr]) * (
                                                 design1[row_nr, factor_nr] + 1) / 2

    return design1, design2


class REGRESSION:
    """
    Regression returns the b0 and b1 coefficients of single or multiple linear regressions.
    It also serves to streamline implementation of stepwise factor minimisation.
    """

    def __init__(self, DataBall_object, x_list: list, y: str):
        self.DataBall_object = DataBall_object
        self.data = DataBall_object.dataset()
        self.x_list = x_list
        self.y = y
        pass

    def linear_reg(self, *args):
        try:
            x_list = args[0]
            ANOVA_object = ANOVA(self.DataBall_object, x_list, self.y)
        except:
            ANOVA_object = ANOVA(self.DataBall_object, self.x_list, self.y)

        if len(self.x_list) == 1:
            self.rc = [ANOVA_object.regression()]  # contains slope and intercept coefficients, respectively

        if len(self.x_list) > 1:
            self.rc = ANOVA_object.coefficient_matrix()  # contains intercept and slope coefficients.

        MODEL_ANALYSIS_object = MODEL_ANALYSIS(ANOVA_object)
        MODEL_ANALYSIS_object.RUN()

        return MODEL_ANALYSIS_object.data

    def AIC_criterion(self, *args):
        try:
            y_pred = args[0]["y_predicted"].values
            print("shape: ", y_pred.shape, "type: ", type(y_pred))
        except:
            y_pred = self.linear_reg()["y_predicted"].values

        try:
            x_list = args[1]
        except:
            x_list = self.x_list

        y = self.data[self.y]

        print(x_list)
        if len(x_list) == 0:
            return dmba.AIC_score(y, [np.average(y.values)] * y.shape[0], df=1)
        return dmba.AIC_score(y, y_pred, df=len(x_list))

    def BIC_criterion(self):
        y = self.data[self.y]
        y_pred = self.linear_reg()["y_predicted"]

        if len(self.x_list) == 0:
            return dmba.BIC_score(y, np.average(y.values) * y.shape[0], df=1)
        return dmba.BIC_score(y, y_pred, df=len(self.x_list))

    def stepwise_minimization(self):
        best_model, best_variables = dmba.stepwise_selection(self.x_list, self.linear_reg, self.AIC_criterion,
                                                             direction="bacwkard", verbose=True)
        return best_model, best_variables


class ANOVA:

    def __init__(self, DataBall_object, x_list: list, y):
        """
        Does ANOVA testing for single and multiple factor models.
        x factor(s) must be passed as a list (e.g. ["cycles","denaturation time", ... ]
        y must be passed a string
        :param DataBall_object:
        :param x_list:
        :param y:
        """

        self.data = DataBall_object.dataset()  # MUST BE A PANDAS DATA FRAME. has the DOE design + it's results.
        self.factor_DF = DataBall_object.factor_DF  # MUST BE A PANDAS DATA FRAME. has all the info on the input factors.

        if type(x_list) == list and len(x_list) > 1:
            self.input_data = self.data[x_list]
            self.factor_label = x_list
        if type(x_list) == list and len(x_list) == 1:
            self.input_data = self.data[x_list]
            self.factor_label = x_list[0]

        self.input_data.insert(0, 'identity', pd.Series([1 for n in range(self.input_data.shape[0])]))
        self.output_data = DataBall_object.dataset()[y]
        self.result_name = y

    # (SIMPLE) LINEAR REGRESSION ---------------------------------------------------
    def SSxy(self):
        # SSfactor = SUM(Xi*Yi) - n Xavg*Yavg

        part1 = sum(np.prod([self.data[self.factor_label], self.data[self.result_name]], axis=0))
        part2 = self.data.shape[0] * np.average(self.data[self.factor_label].values) * np.average(
            self.data[self.result_name].values)
        SSFactor = part1 - part2
        return SSFactor

    def SS(self, column):
        # SS = SUM(Xi^2) - n * Xavg^2
        part1 = sum([self.data.loc[row, column] ** 2 for row in range(self.data.shape[0])])
        part2 = self.data.shape[0] * np.average(self.data.loc[:, column]) ** 2
        SS = part1 - part2
        return SS

    def slope(self):
        return self.SSxy() / self.SS(self.factor_label)

    def SSregression(self):
        return self.slope() * self.SSxy()

    def SStotal(self):
        # SUM(Yi^2) - SUM(Yi)^2/n
        part1 = sum(self.data.loc[:, self.result_name] ** 2)
        part2 = (sum(self.data.loc[:, self.result_name]) ** 2) / self.data.shape[0]
        SStotal = part1 - part2
        return SStotal

    def SSerror(self, type="simple"):
        if type == "simple":
            return self.SStotal() - self.SSregression()
        if type == "matrix":
            return self.SStotal_matrix() - self.SSregression_matrix()

    def degreesFreedom(self, type):
        if type == "regression":
            return 1
        elif type == "error":
            return self.data.shape[0] - 2
        elif type == "total":
            return self.data.shape[0] - 1

    def MSquare(self, type):
        # type can be "regression" or "error"
        if type == "regression":
            return self.SSregression() / self.degreesFreedom(type=type)
        elif type == "error":
            return self.SSerror() / self.degreesFreedom(type=type)

    def Ftest(self):
        return self.MSquare("regression") / self.MSquare("error")

    def ANOVA_test(self):
        header = ["Sum of Squares", "df", "Mean Square", "F ratio"]
        rows = ["Regression", "Error", "Total"]
        data = []

        for i in range(len(rows)):
            if i == 0:
                data.append(
                    [self.SSregression(), self.degreesFreedom(rows[i].lower()), self.MSquare(rows[i].lower()),
                     self.Ftest()])
                continue
            elif i == 1:
                data.append([self.SSerror(), self.degreesFreedom(rows[i].lower()),
                             self.MSquare(rows[i].lower())])
                continue
            elif i == 2:
                data.append([self.SStotal(), self.degreesFreedom(rows[i].lower())])
                continue
        return pd.DataFrame(data, index=rows, columns=header)

    def regression(self):
        slope = self.SSxy() / self.SS(self.factor_label)
        intercept = np.average(self.data.loc[:, self.result_name]) - slope * np.average(
            self.data.loc[:, self.factor_label])
        return slope, intercept

    def t_test(self, slope_coefficient, intercept_coefficient):
        t0_intercept = intercept_coefficient / np.sqrt(self.MSquare("error") * self.data.shape[0] + (
                np.average(self.data.loc[:, self.factor_label]) / self.SS(self.factor_label)))
        t0_slope = slope_coefficient / np.sqrt(self.MSquare("error") / self.SS(self.factor_label))

        intercept_p_value = 2 * np.e ** (-t0_intercept ** 2 / 2) / np.sqrt(2 * np.pi)
        slope_p_value = 2 * np.e ** (-t0_slope ** 2 / 2) / np.sqrt(2 * np.pi)

        return intercept_p_value, slope_p_value

    # MULTIPLE LINEAR REGRESSION ---------------------------------------------------
    def coefficient_matrix(self):
        X = self.input_data.values  # Changing which x_list are passed via self.input_data will return different slopes.

        Y = self.output_data.values
        XtX = X.T @ X
        b1 = np.linalg.pinv(XtX) @ X.T @ Y

        # b1 = np.linalg.pinv(self.input_data.transpose() * self.input_data) * self.input_data.transpose()* self.out_data.iloc[0]

        return b1

    def SStotal_matrix(self):
        Y = self.output_data.values
        YtY = Y.T @ Y
        nY2 = len(Y) * np.average(Y) ** 2

        return YtY - nY2

    def SSregression_matrix(self):
        X = self.input_data.values
        Y = self.output_data.values
        b1 = self.coefficient_matrix()

        B1XtY = b1.T @ X.T @ Y
        nY2 = len(Y) * np.average(Y) ** 2

        return B1XtY - nY2

    def SSerror_matrix(self):
        return self.SStotal_matrix() - self.SSregression_matrix()

    def degreesFreedom_matrix(self, type):
        k = self.input_data.shape[1]
        n = self.input_data.shape[0]
        if type == "regression":
            return k - 1
        if type == "error":
            return n - k  # k is 12 + 1. 12 factors + 1 identity column. the - 1 here is omitted.
        if type == "total":
            return n - 1

    def MSquare_matrix(self):
        SSregression = self.SSregression_matrix()
        SStotal = self.SStotal_matrix()
        SSerror = SStotal - SSregression
        df_r = self.degreesFreedom_matrix("regression")
        df_e = self.degreesFreedom_matrix("error")
        return SSregression / df_r, SSerror / df_e

    def Ftest_matrix(self):
        a, b = self.MSquare_matrix()
        return a / b

    def Rsquared_matrix(self):
        Rsquared = self.SSregression_matrix() / self.SStotal_matrix()
        Rsquared_adjusted = 1 - (self.degreesFreedom_matrix("total")) / (self.degreesFreedom_matrix("error")) * (
                1 - Rsquared)
        return [Rsquared, Rsquared_adjusted]

    def MANOVA_test(self):
        header = ["Sum of Squares", "df", "Mean Square", "F ratio", "R2", "R2_adjusted"]
        rows = ["Regression", "Error", "Total"]
        data = []

        for i in range(len(rows)):
            if i == 0:
                data.append(
                    [self.SSregression_matrix(), self.degreesFreedom_matrix("regression"), self.MSquare_matrix()[0],
                     self.Ftest_matrix(), self.Rsquared_matrix()[0], self.Rsquared_matrix()[1]])
                continue
            elif i == 1:
                data.append([self.SSerror_matrix(), self.degreesFreedom_matrix("error"),
                             self.MSquare_matrix()[1]])
                continue
            elif i == 2:
                data.append([self.SStotal_matrix(), self.degreesFreedom_matrix("total")])
                continue
        return pd.DataFrame(data, index=rows, columns=header)


class MODEL_ANALYSIS():
    def __init__(self, ANOVA_object):
        self.rc = ANOVA_object.coefficient_matrix()
        self.data = ANOVA_object.data
        self.input_data = ANOVA_object.input_data.copy()
        self.out_data = ANOVA_object.output_data.copy()
        self.y = ANOVA_object.result_name

    def hat_matrix(self):
        data = self.data.values
        return data @ np.linalg.pinv(data.T @ data) @ data.T

    def y_pred(self, append=True):
        row_count = self.data.shape[0]
        coefficients = self.rc
        y_pred = []

        for i in range(row_count):
            row_values = self.input_data.iloc[i, :]
            y_pred.append(sum(row_values * coefficients))

        if append:
            self.data["y_predicted"] = y_pred
        elif not append:
            return y_pred

    def resid(self, append=True):
        if append:
            self.data["residual"] = self.data["y_predicted"] - self.out_data
        elif not append:
            return [self.data.loc["y_predicted"] - self.out_data]

    def MSE(self, column: str):
        n = self.data.shape[0]
        SUMe2 = sum(self.data[column] ** 2)
        return SUMe2 / n

    def std_resid(self, append=True):
        MSEresiduals = self.MSE("residual")
        hat_matrix = np.array(self.hat_matrix())
        standardized_residuals = []
        row_nr = self.data.shape[0]

        # ith row standard_residual = (ith row residual) / sqrt(Mean Squared error of residuals * (1 - diagonal of hat matrix))
        for i in range(row_nr):
            ith_std_resid = self.data.loc[i, "residual"] / np.sqrt(MSEresiduals * (1 - hat_matrix[i, i]))
            standardized_residuals.append(ith_std_resid)

        if append:
            self.data["std_residual"] = standardized_residuals
        if not append:
            return standardized_residuals

    def RUN(self):
        self.y_pred()
        self.resid()
        # self.std_resid()

def plot_stacking(function):
    """
     Allows a given plotting function to overlay other data. data is passed in via dataset=<DataBall.dataset()>
     the function must have **kwargs.
     code that sets the working dataset must be wrapped into accepting data from kwargs["dataset"]
     the overall plot function must return fg and ax.

     :param function:
     :return:
     """
    def wrapper(self,*args, **kwargs):
        fg, ax = function(self, *args)
        for i in kwargs["add"]:
            fg, ax = function(self, *args, fg=fg, ax=ax, dataset=i)
        plt.show()
    return wrapper


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

        self.factor_names = HTMLDataExtractor(self.continuousFactors,
                                              "name")  # extracts names from the website but these are different from the source documentation
        self.factor_attributes = ["min", "max", "value", "type"]
        self.factor_min = HTMLDataExtractor(self.continuousFactors, "min")
        self.factor_max = HTMLDataExtractor(self.continuousFactors, "max")
        self.factor_setvalues = HTMLDataExtractor(self.continuousFactors, "value")
        self.factor_datatype = HTMLDataExtractor(self.continuousFactors, "type")
        self.factor_polymerase = HTMLDataExtractor(self.continuousFactors, "value")
        self.factor_class_integers = {"id_cycles": 0}
        self.datapackage = [self.factor_min, self.factor_max, self.factor_setvalues, self.factor_datatype]

        self.factor_DF = dataframe_generate(self.datapackage, self.factor_attributes, self.factor_names)
        self.factor_names = [factors for factors in self.factor_DF.columns]  # updates factor names with the true ones.
        self.factor_count = len(self.factor_names)
        self.factor_DF_cache = []  # Structure: [ [index:int,"DOE_<version> - <DOE_design> factor_DF <version>", NxM matrix:list of lists] ]

        self.DOE_cache = []  # Structure: [ [index:int,"DOE_<version> - <DOE_design>", NxM matrix:list of lists] ]
        self.DOE_version = 0
        self.DOE_active_pointer = 0

        self.output_labels = ["DNA yield (ng/uL)", "DNA Purity (%)", "x amplification", "duration (s)"]
        self.output_active_pointer = 0

        self._TEMP_DATA = False  # gets set to true when modifying data. Makes the self.plot_...() functions use self.DOE_active, instead of cached DOE data.
        self._LAZY_OUTPUT_SELECTION = False  # gets set to true when modifying data. remembers the last output column

        self._LOCAL_STORAGE = {}  # cache for names and versions of imported DOEs / Structure {"name":str : version:int}

    # Decorators ------------------------------------------

    # Saving and Importing functions -------------------------------------------------------------------
    def _folder(self):
        year, month, day, *_ = time.localtime()
        foldername = f"{year}_{month}_{day} - Results Folder"
        try:
            os.makedirs(foldername)
            return foldername
        except FileExistsError:
            return foldername

    def _filename(self):
        return self.DOE_cache[self.DOE_active_pointer - 1][1] + ".csv"

    def exportToDirectory(self):
        old_dir = os.getcwd()
        os.chdir(self._folder())

        FILENAME = self._filename()
        dataframe = self.dataset()

        try:
            pd.read_csv(FILENAME)
            with open(FILENAME, "a", newline="") as file:
                writer = csv.writer(file)

                for key, value in enumerate(dataframe.values):
                    writer.writerow(value)

        except FileNotFoundError:
            dataframe.to_csv(FILENAME, index=False)

        os.chdir(old_dir)

    def importFromDirectory(self):
        """
        Grabs the FILE.csv in the local directory and creates a new instance in DataBall.
        the new instance will be labelled self.FILE_import<version>
        <version> is pulled from self._version, where the program checks for duplicate imports and updates the version.

        :return: [.csv file file_name, .csv data]
        """

        old_directory = os.getcwd()

        def openFile():
            while True:
                for keys, files in enumerate(os.listdir()):
                    print(f"keys: {keys}, files: {files}")

                selection = input("select the index of a folder or .csv file to open. (type .. to go back in folders)")

                try:  # goal here is to generate the filename and csv objects.
                    filename = os.listdir()[int(selection)]
                    CSV = pd.read_csv(filename)
                    print("{} opened successfully".format(filename))
                    break
                except:
                    pass

                try:
                    os.chdir(os.listdir()[int(selection)])
                except:
                    pass

                if selection == "..":
                    os.chdir(selection)
            return CSV, filename

        CSV, filename = openFile()  # seeks and opens .CSV file as pandas.DataFrame. Current directory has been changed.

        def fileVersion(newname):
            version = str(next(self._version(newname)))
            return version

        def nameFile():
            name = os.path.splitext(filename)[0]
            suffix = "_import"
            newname = name + suffix
            return name, newname

        name, newname = nameFile()  # creates the new file_name for the imported file
        version = fileVersion(newname)  # versions the imported file.
        os.chdir(old_directory)  # resets the directory

        def dataStructure(df_object):
            col = df_object.columns

            print("\n \n  --- file structuring --- ")
            print("--- BASIC FILE STRUCTURE:  --- \n")
            print("1ST ROW (cell A1): columns labels // COLUMN ORDER: inputs - outputs - analysis")
            print("OTHER ROWS: column-specific values // COLUMN ORDER: inputs - outputs - analysis \n")
            test_factor_count = int(input("How many test factor columns (or inputs) does the file have?"))
            self.factor_count = test_factor_count
            self.factor_names = list(col[:self.factor_count])

            output_factor_count = int(input("How many result columns (or outputs) does the file have?"))
            self.output_count = output_factor_count
            self.output_labels = list(col[self.factor_count: self.factor_count + self.output_count])

            analysis_column_count = len(col) - test_factor_count - output_factor_count
            analysis_columns = list(col[self.factor_count + self.output_count:])

            columns = self.factor_names + self.output_labels + analysis_columns
            return columns

        columns = dataStructure(CSV)  # sets the structure for what is an input, an output and analysis data.

        def datatableGeneration(name, version, data, column_labels):
            setattr(self, name + version, pd.DataFrame(update_DOEmatrix_datatypes_int64(
                data.values, **self.factor_class_integers), columns=column_labels))

        datatableGeneration(newname, version, CSV,
                            columns)  # instances a variable in .self with the specified file_name and data.
        print(f"{name} was imported as self.{newname}{version}")

        def datatableSelection():
            self.DOE_active = getattr(self, newname + version)
            self.DOE_active = self.DOE_active.astype({"cycles": int})
            self.DOE_active = self.DOE_active.astype({"cycles": "float64"})

        datatableSelection()  # sets active DOE datatable as a copy of the instanced variable.

        def cacheUpdate():
            self.DOE_active_pointer = len(self.DOE_cache)  # sends the pointer to the top of the cache.
            self.DOE_active_pointer += 1  # updates the pointer.
            self.DOE_cache.append([version, newname + version, getattr(self, newname + version)])  # adds the import

        cacheUpdate()  # puts a copy of the newly instanced variable into a cache for future calls.

        factor_DF_filename = "factor_DF_" + filename

        def factorInfoGeneration():
            columnIndex = self.factor_count
            columnLabels = self.dataset().columns[:columnIndex]
            rowLabels = ["min", "max", "value"]

            # structure: [[min1,max1,0],[min2,max2,0],...]
            values = [[min(self.dataset().iloc[:, index]), max(self.dataset().iloc[:, index]), 0] for index in
                      range(columnIndex)]

            # structure: [[min1,min2,...],[max1,max2,...],[0,0,0...]]
            values_T = zip(*values)

            dataframe = pd.DataFrame(values_T, index=rowLabels, columns=columnLabels)

            setattr(self, factor_DF_filename, dataframe)

        factorInfoGeneration()  # creates the factor_DF dataframe file and makes it an instance in .self

        def factorDFSelection(file_name):
            self.factor_DF = getattr(self, file_name)

        factorDFSelection(factor_DF_filename)  # sets self.factor_DF = dataframe

        def factorDFCacheUpdate():
            factor_DF_version = fileVersion(factor_DF_filename)  # versions the factor_DF
            self.factor_DF_cache.append([factor_DF_version, factor_DF_filename, self.factor_DF])

        factorDFCacheUpdate()

    # Pointer-guided information retrieval functions -------------------------------------------------------------------
    def dataset(self):  # returns what the active pointer is looking at
        if self._TEMP_DATA: return self.DOE_active
        if not self._TEMP_DATA: return getattr(self, self.DOE_cache[self.DOE_active_pointer - 1][1])

    def dataset_DOEmatrix(self):
        if self._TEMP_DATA: return self.DOE_active.iloc[:, :self.factor_count]
        if not self._TEMP_DATA: return getattr(self, self.DOE_cache[self.DOE_active_pointer - 1][1]).iloc[:,
                                       :self.factor_count]

    def dataset_results(self):
        if self._TEMP_DATA: return self.DOE_active.iloc[:, self.factor_count:]
        if not self._TEMP_DATA: return getattr(self, self.DOE_cache[self.DOE_active_pointer - 1][1]).iloc[:,
                                       self.factor_count:]

    def dataset_reset(self):
        self.ANALYTICS_mode()
        self.DOE_active = getattr(self, self.DOE_cache[self.DOE_active_pointer - 1][1])

    def dataset_results_pointer(self):
        if self._LAZY_OUTPUT_SELECTION: return self.DOE_active.iloc[:, self.output_active_pointer]
        if not self._LAZY_OUTPUT_SELECTION: return getattr(self, self.DOE_cache[self.DOE_active_pointer - 1][1]).iloc[:,
                                                   self.output_active_pointer]

    #  Generators for state changes --------------------------------------------------------------------------------
    def _version(self, filename: str):
        """
        Returns the next version number for an imported DOE design.

        filename will be in the format of <DOE_design><version>_import (e.g. sobol1_data_import)

        :param filename:
        :return: sobol1_data_import1 if unique, or import2,3,4,5... if duplicate
        """

        VERSION = 1

        while True:
            if filename not in self._LOCAL_STORAGE.keys():
                self._LOCAL_STORAGE[filename] = VERSION
            elif filename in self._LOCAL_STORAGE.keys():
                self._LOCAL_STORAGE[filename] += 1
            yield self._LOCAL_STORAGE[filename]

    def _objectStates(self, name, states, return_cache=0):
        GENERATOR_CACHE = "_GENERATOR_CACHE"
        if return_cache == 1:
            return getattr(self, GENERATOR_CACHE)

        try:
            getattr(self, GENERATOR_CACHE)
        except AttributeError:
            setattr(self, GENERATOR_CACHE, {})

        try:
            getattr(self, GENERATOR_CACHE)[name]
        except KeyError:
            getattr(self, GENERATOR_CACHE)[name] = itertools.cycle(states)

    # DataBall mode switches -------------------------------------------------------------------------------
    def ANALYTICS_mode(self):
        name = "ANALYTICS"
        states = ["ON", "OFF"]
        GENERATOR_CACHE = "_GENERATOR_CACHE"

        self._objectStates(name, states)
        state = next(getattr(self, GENERATOR_CACHE)[name])

        if state == states[0]:  # i.e. "ON"
            self._TEMP_DATA = True
            self._LAZY_OUTPUT_SELECTION = True
            print(f"analytics mode:{state}")
            return state
        elif state == states[1]:  # i.e. "OFF"
            self._TEMP_DATA = False
            self._LAZY_OUTPUT_SELECTION = False
            print(f"analytics mode:{state}")
            return state

    # Input mapping function (Maps variables to pyDOE's DoE function format) -----------------------------------
    def _function_mapping(self, function):
        function_map = \
            {sobol: [self.factor_DF.loc["min"], self.factor_DF.loc["max"], int(input("How many runs for the sobol?")),
                     self.factor_DF.shape[1]],
             ff2n: [self.factor_count, self.factor_DF.loc["min"], self.factor_DF.loc["max"]],
             fractionalFactorial: [self.factor_count, 4, self.factor_DF.loc["min"], self.factor_DF.loc["max"]]
             }
        # factors, resolution, min, max
        return function_map[function]

    # DOE design generation, caching, updating & retrieval and simulation running functions ------------------------
    def DOE_import(self, design, self_data=True, **kwargs):
        """
        Takes in a function via the design variable.
        Passes the function into function_mapping() to collect the arguments to pass into itself.
        Creates self."design<version>" variable to hold the updated DoE matrix, in int64 format.
        Creates self."design<version>.1" variable to hold the boilerplate DoE matrix.
        Caches the design and primes the DoE into the system (via pointer)

        if self_data=False: (When you're importing data)


        """

        DOE_version = str(next(self._version(design.__name__)))

        if self_data:
            data_package = self._function_mapping(design)
            setattr(self, design.__name__ + DOE_version,
                    update_DOEmatrix_datatypes_int64(design(data_package)[0], **self.factor_class_integers))
            setattr(self, design.__name__ + DOE_version + "a", design(data_package)[1])

        print(self.__class__.__name__, ":", design.__name__ + DOE_version, "generated")
        print(self.__class__.__name__, ":", design.__name__ + DOE_version + "a",
              "generated ---> blank matrix")

        self.DOE_active_pointer = len(self.DOE_cache)  # sends the pointer to the top of the cache.
        self.DOE_active_pointer += 1  # adds 1 to mark the entry of a new DOE design.
        self.DOE_cache.append(
            [DOE_version, design.__name__ + DOE_version,
             pd.DataFrame(getattr(self, design.__name__ + DOE_version), columns=self.factor_names)])
        self.DOE_active = pd.DataFrame(getattr(self, design.__name__ + DOE_version), columns=self.factor_names)

        self.DOE_active = self.DOE_active.astype({"cycles": int})
        self.DOE_active = self.DOE_active.astype({"cycles": "Int64"})

        if not self_data:
            pass

    def DOE_current_design(self, change=""):
        """
        Shows all  DOE designs currently stored in the cache of this DataBall class object.
        change: int allows the user to select which DOE design and it's underlying data gets primed for the simulation.

        :change: if "" > returns all DOE designs by their order of generation
        :change: if int > moves self.DOE_active_pointer to stated DOE design, alongside all underlying data.

        :return: see above comments for :change:
        """
        if change == "":
            for row_nr in range(len(self.DOE_cache)):
                print("DOE nr.: {0} // name: {1}".format(row_nr + 1, self.DOE_cache[row_nr][1]))
            print("Active DOE: >>{}<<".format(self.DOE_active_pointer))
            return

        try:
            change_int = int(change)
        except ValueError:
            print("Invalid input: expected an integer or empty string")
            return

        try:
            self.DOE_active = self.DOE_cache[change_int - 1][2]  # DOE cache has the matrix in index 2.
            self.DOE_active = self.DOE_active.astype({"cycles": int})
            self.DOE_active = self.DOE_active.astype({"cycles": "Int64"})
            self.DOE_active_pointer = change_int

            self.factor_DF = self.factor_DF_cache[change_int - 1][2]  # factor_df cache has the matrix in index 2.

        except IndexError:
            print("That index is out of bounds")

    def DOE_update(self, change=False):

        transposed_DF = self.factor_DF.transpose()

        if not change:
            print(self.factor_DF.transpose())
        if change:
            print("Factor range updates:")
            for index in self.factor_DF.transpose().index:
                for column in self.factor_DF.transpose().columns[[0, 1]]:
                    try:
                        transposed_DF.loc[index, column] = int(input(f"{index} {column}:"))
                    except:
                        transposed_DF.loc[index, column] = transposed_DF.loc[index, 2]

        self.factor_DF = transposed_DF.transpose()

    def RUN(self, hard_limit=None):
        """
        Runs Ben's PCR simulator.
        Runs the DOE design assigned by the pointer.

        :param recent_design:
        :param args:
        :return:
        """
        while True:
            if self.ANALYTICS_mode() == "OFF":
                break

        results = []
        DOE_matrix = self.DOE_active
        try:
            if type(int(hard_limit)) == int:
                DOE_matrix = DOE_matrix.iloc[:int(hard_limit), :]  # truncates the DoE design
        except:
            pass
        for rows in DOE_matrix.values:
            test_run = dataframe_to_pcr_format(rows)
            results.append(pcrsim.demo(test_run))  # This takes forever and it needs to be sorted
            print("{0} out of {1}".format(len(results), len(DOE_matrix)))

        table_DF = self._DOE_extract_data(results)

        setattr(self, self.DOE_cache[self.DOE_active_pointer - 1][1], table_DF)
        self.DOE_active = self.dataset()

    # data wrangling functions -----------------------------------------------------
    def _DOE_extract_data(self, data):

        interim_df1 = pd.DataFrame(data)
        interim_df2 = pd.DataFrame(interim_df1.values, columns=self.output_labels)
        final_df = self.DOE_active.join(interim_df2)

        return final_df

    def data_topvalues(self, n, column_selection, sort_ascending=False, object=False):

        data = self.data_sort(column_selection, sort_ascending=sort_ascending)

        if self._TEMP_DATA == True:
            self.DOE_active = self.DOE_active.iloc[:n, :]

            if object == True:
                return self.DOE_active

        elif self._TEMP_DATA == False or object == True:
            return data.iloc[:n, :]

    def data_average_ranking(self):
        pointer_state = self._TEMP_DATA

        rankings = []

        data = self.dataset()
        column_index = [12, 13, 14, 15]
        _temp_attribute_list = [f"df{nr}" for nr in range(len(self.output_labels))]

        version = 0

        self._TEMP_DATA = False
        for i in column_index:
            self.output_active_pointer = i
            sorted_column = self.data_sort(sort_ascending=False, object=True).iloc[:, i]
            setattr(self, f"df{version}", sorted_column)
            version += 1

        for df in _temp_attribute_list:
            _df = getattr(self, df)
            for i in range(len(data)):
                _df.iloc[i] = i

        for row in range(len(data)):
            result = []
            for df in _temp_attribute_list:
                result.append(getattr(self, df).loc[row])
            rankings.append(result)

        rankings = pd.DataFrame(rankings, columns=self.output_labels)
        self.DOE_active["Average Ranking"] = rankings.mean(axis=1)

        self._TEMP_DATA = pointer_state

    def data_sort(self, column_selection: str, sort_ascending=True, object=False):

        try:
            self.output_active_pointer = self.dataset().columns.get_loc(column_selection)
        except:
            pass

        if self._LAZY_OUTPUT_SELECTION:
            selection = self.output_active_pointer

        else:
            for i in range(len(self.output_labels)):
                print("{}. - {}".format(i, self.output_labels[i]))
            while True:
                try:
                    selection = int(input("Select an output to sort by")) + self.DOE_active.shape[1]
                    break
                except:
                    print("selection needs to be an integer")

        if self._TEMP_DATA:
            self.DOE_active.sort_values(self.DOE_active.columns[selection], axis=0, inplace=True,
                                        ascending=sort_ascending)
            if object == True:
                return self.DOE_active

        elif self._TEMP_DATA == False or object == True:
            return self.DOE_active.sort_values(self.DOE_active.columns[selection], axis=0, ascending=sort_ascending)

    # data plotting functions -----------------------------------------------------
    @plot_stacking
    def plot_byoutput(self,**kwargs):

        try:
            data = kwargs["dataset"]
        except:
            data = self.dataset()

        try:
            fg, axs = [kwargs["fg"], kwargs["ax"]]
        except:
            fg, axs = plt.subplots(nrows=2, ncols=2, figsize=(5.5, 3.5), layout="constrained")

        x_axis = [str(i) for i in range(data.shape[0])]

        # self.data() returns the most recent a.RUN() results (DataFrame object)
        # self.data().count() returns a list with the count of non-NaN values in each column
        # self.data().count()[self.output_labels[0]] returns the count of non-NaN values for column = "DNA yield (ng/uL)"
        # x_axis = a list with the rows which got populated in RUN()

        grid = [position for position in itertools.product([0, 1], repeat=2)]
        label = self.output_labels

        plot_layout = zip(grid, label)

        for plot in plot_layout:
            axs[plot[0]].scatter(x_axis, data[plot[1]], marker=".")
            axs[plot[0]].set_title(plot[1])

        return fg, axs

    def plot_byfactor(self, *regression_data):

        data = self.dataset()
        try:
            best_variables = regression_data[0]
            coefficient_vector = regression_data[1]
            coefficient_vector2 = coefficient_vector.reshape(1, len(coefficient_vector))  # weird bug?
            averages = [0] + [np.average(data[factor]) for factor in best_variables]
            d = [coefficient_vector2.flatten(), averages]

            coef_DF = pd.DataFrame(d, columns=["identity"] + best_variables)

            _temp_DF = coef_DF

        except:
            print("regression table not made")

        counter = 0

        fg, axs = plt.subplots(nrows=4, ncols=3, figsize=(10, 10), layout="constrained")

        DOE_matrix, results = data.iloc[:, :self.factor_count], data.iloc[:, self.factor_count:]
        # 0 - self.factor_count = DOE_matrix
        # self.factor_count - 16 = results

        y_axis = results.iloc[:, 0].values  # yield column

        for row in range(len(axs)):
            for column in range(len(axs[row])):
                x_axis = DOE_matrix.iloc[:, counter].values
                axs[row, column].scatter(x_axis, y_axis)
                axs[row, column].set_title(self.factor_names[counter])
                factor = self.factor_names[counter]
                try:
                    if factor in coef_DF:
                        # plots x_axis against x_axis * slope + intercept + average response of other factors

                        _temp_DF = coef_DF.drop(columns=[factor, "identity"])

                        _avg = (_temp_DF.loc[0] * _temp_DF.loc[1]).sum()  # Sum of average results of all other factors

                        y_pred = x_axis * coef_DF.loc[0, factor] + coef_DF.loc[0, "identity"] + _avg

                        axs[row, column].plot(x_axis, y_pred, color="0")
                        _temp_DF = coef_DF
                except:
                    pass
                counter += 1
        plt.show()

    @plot_stacking
    def plot_distributions(self, column, **kwargs):
        """
        Allows a given plotting function to overlay other data. data is passed in via dataset=<DataBall.dataset()>
        the function must have **kwargs.
        code that sets the working dataset must be wrapped into accepting data from kwargs["dataset"]
        the overall plot function must return fg and ax.

        :param function:
        :return:
        """
        try:
            data = kwargs["dataset"]
            print("imported new dataset")
            print(data)
        except:
            data = self.dataset()
            print("using current dataset")

        if self._LAZY_OUTPUT_SELECTION:
            selection = self.output_active_pointer
        else:
            selection = data.columns.get_loc(column)

        cumulative_frequency = [(row - 0.5) / len(data) for row in range(data.shape[0])]
        normal_frequency = [st.norm.ppf(p) for p in cumulative_frequency]

        try:
            fg, axs = [kwargs["fg"], kwargs["ax"]]
        except:
            fg, axs = plt.subplots(nrows=1, ncols=4, figsize=(11, 3.5), layout="constrained")

        axs[0].scatter(data.iloc[:, selection].values, normal_frequency)
        axs[0].set_title("Normal Probability plot")
        axs[0].set_xlabel(column)
        axs[0].set_ylabel("Z-score")

        axs[1].hist(data.iloc[:, selection].values)
        axs[1].set_title("Data distribution")
        axs[1].set_xlabel(column)
        axs[1].set_ylabel("counts")
        axs[1].axvline(x=0, color="red", linewidth="1", linestyle="dashed")

        def skip_row(n):
            try:
                return data.loc[n][column]
            except:
                print("skipping over removed rows")
                pass

        factor_by_original_row_index = [skip_row(i) for i in range(data.shape[0])]

        row_index = [i for i in range(data.shape[0])]
        axs[2].scatter(row_index, factor_by_original_row_index)
        axs[2].set_title("{} by Row".format(column))
        axs[2].set_xlabel("Row nr.")
        axs[2].set_ylabel(column)
        axs[2].axhline(y=0, color="red", linewidth="1", linestyle="dashed")

        axs[3].scatter(data[column], data["y_predicted"])
        axs[3].set_title("y_predicted by {}".format(column))
        axs[3].set_xlabel(column)
        axs[3].set_ylabel("y_predicted")
        axs[3].axvline(x=0, color="red", linewidth="1", linestyle="dashed")

        return fg, axs

    @plot_stacking
    def plot_withingroups(self, statistic=0,**kwargs):
        """
        Plots a statistic (e.g.0 for average,1 for std. dev) for each treatment level of a given test factor (i.e. a DOE mean, or DOE std.dev plot)
        :param statistic: "average", "std. dev."
        :return:
        """
        if statistic == 0:
            stat = np.average
        elif statistic == 1:
            stat = np.std

        try:
            dataset = kwargs["dataset"]
        except:
            dataset = self.dataset()

        factors = dataset.columns[:self.factor_count]
        treatments = dataset.iloc[:, :self.factor_count]
        results = "DNA yield (ng/uL)"

        # treatment_levels = [treatments.value_counts(factor).shape[0] for factor in factors]
        # level_counts = [treatments.value_counts(factor).values for factor in factors]

        # generator keys are factors. for each key, there is a list of length 2. index 0 is the nr. of treatments. index 1 is the nr of rows in each treatment
        generator = {factor: [treatments.value_counts(factor).shape[0], treatments.value_counts(factor).values] for
                     factor in factors}

        def average_by_level(generator):
            averages = {}

            for factor in generator.keys():

                dataset.sort_values(factor, axis=0, ascending=True, inplace=True)
                factor_avg = []
                start_index = 0

                for treatment in range(generator[factor][
                                           0]):  # treatments is an integer, generator[factor][0] is how many treatments there are.

                    subset_rows = generator[factor][1][treatment]

                    value = stat(dataset[results][start_index:start_index + subset_rows])  # average of subset of rows.
                    factor_avg.append(value)
                    start_index += subset_rows  # shifts the start point of the rows.

                averages[factor] = factor_avg

            return averages

        factor_averages = average_by_level(generator)  # {factor1: [50], factor2: [23,45] ....}

        try:
            fg, ax = [kwargs["fg"],kwargs["ax"]]
        except:
            fg, ax = plt.subplots(figsize=(7, 4.5), layout="constrained")

        def bars_by_group(averages):
            position = 0.5
            for key in averages:
                increment_value = 1 / (1 + len(averages[key]))
                x_coordinates = [[position + i * increment_value] for i in range(1, 1 + len(averages[key]))]

                ax.plot(x_coordinates, averages[key], mec="k", marker="x")

                position += 1.5

        bars_by_group(factor_averages)
        x = [1 + 1.5 * i for i in range(len(factors))]
        ax.set_xticks(x, labels=factors, rotation=90, rotation_mode="default")

        ax.axhline(stat(dataset[results]), color="red", linewidth=1, linestyle="dashed")

        if statistic == 0:
            ax.set_title("DOE plot - Average per treatment level ({})".format(results))
            ax.set_ylabel("Average")
        elif statistic == 1:
            ax.set_title("DOE plot - Standard deviation per treatment level ({})".format(results))
            ax.set_ylabel("standard deviation")

        for i in x:
            ax.axvline(i, color="gray", linewidth=0.1)

        return fg,ax


# Decorators ---------------------------------------------------


"""
def plot_stacking(function):
    def wrapper(*args,**kwargs):
        fg,ax = function(*args)
        for i in kwargs["add"]:
            fg,ax = function(i[0],i[1],fg=fg,ax=ax)
        plt.show()
    return wrapper

@plot_stacking
def test_plot(a,b,**kwargs):
    x = np.linspace(1,a,num=10)
    y = np.geomspace(1,b,num=10)
    try:
        fg, ax = [kwargs["fg"],kwargs["ax"]]
    except:
        fg, ax = plt.subplots()
        
    ax.scatter(x,y)
        
    return fg, ax

test_plot(4,400,add=[[2,50],[4,100]])

"""

## script execution ----------------------------------

a = DataBall()  # creates the DataBall object
# a.DOE_import(sobol)  # creates an DOE design
# a.RUN(hard_limit=3)  # runs the DOE design
