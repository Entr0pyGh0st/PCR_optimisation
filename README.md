

Experimental designs are methods of querying the behaviour of complex in silico/in vivo systems in an extensive and data-efficient manner.

They address a technical execution problem that arises when a scientist/engineer seeks to investigate and model a black box system with large sets of input parameters empirically. Generally speaking, both the inputs and output are usually continous and the resulting approximation model is usually a simple combination of multiple linear or quadratic models. 

Whilst individual input parameter testing is feasible, accepteable and of straight forward analysis and experimental control, it's cumbersome nature of iterating through each individual factor is both resource ineficient and "low resolution" (i.e. doesn't explain how inputs can interact with eachother to impact the output, only how each input affects the output exclusively).   

A subset of experimental designs addresses this by providing bespoke testing schematics that leverage controlled input aliasing (i.e. having two factors changing simultaneously in a test run) to 1) to test more inputs per test run and 2) investigate for input interactions (n-factor interactions), thus allowing for a more complete description of black box behaviour compared straightforward indididual factor testing. 

Another subset of experimental designs use quasirandom sequences of numbers (low-discrepancy sequences) to generate a testing schematic that covers the experimental ranges of the input factors of a black box evenly. These are technically superior to the latter subset if not constrained by resources or time.

Ben Shirt-Ediss made https://virtual-pcr.ico2s.org/pcr/, an in silico model of a PCR reaction aimed at amplifying a 1kb DNA sequence by changing the parameterization of 12 different inputs. The source code can be found at https://bitbucket.org/ben_s_e/virtual-pcr-notebook/src/main/ .Understandeably, the fact that he coded a model means that one could derive the maximum output metrics (yield, product purity) numerically. However, he and I (as found in this repo), intended that the solution be found via experimental designs, for the sakes of learning about them.

This repo seeks to interface his model with pyDOE, a python module for experimental designs, to provide an interactive approach to experimental design testing without going through the web server. 

# Current highscore: 1.483 mg/mL, 97.6% pure, 8801-fold amplification, 22060 second run 


# Bugs:


# Future features:

more designs:  PB, BB

more plots: interaction plots , Q-Q


# How to use:
Data Generation - Generating a fresh dataset:
```python
instance = DataBall() # creates DataBall class object
instance.DOE_import(sobol) # imports a DOE design, in this case from the function named sobol. Creates an internal attribute named self.<design><version>.
~~~ input request: # requests user input for design-specific variables. for sobol, it will ask for total nr. of runs.
instance.RUN() # runs the newly created attribute ,e.g. self.sobol1,  through the PCR simulator

DOE_update(True) followed by DOE_import(sobol) # 1st function call updates the min,max values. 2nd function call generates an updated sobol design.

```

Data Generation - DOE designs:
```python
Different DOE designs can be readily imported via instance.DOE_import(). To do this, you must pass the function as a keyword into instance.DOE_import().

instance.DOE_import(sobol) - Imports a quasi-random low-discrepancy sobol sequence. number of runs must be specified.
instance.DOE_import(fractionalFactorial2lvl) - Imports a 2 level fractional Factorial design. The resolution of the design must be specified. Note that the function generates the minimum viable design to reach a specified resolution, i.e. uses the minimum amount of main factors to achieve aliasing of the remainder of the factors.
instance.DOE_import(fullFactorial2lvl) - imports a 12-factor full factorial design. The goal is not to run the full design, it's to be able to slice the DOE matrix at your convenience to get full resolution for some factors. slicing is only currently possible by self assignment over instance.dataset().


```


Data Exporting - Saving a dataset:
```python
instance.exportToDirectory() # generates a "year_month_date - Results Folder" in curdir and saves results as sobol1.csv

# sobol1.csv is a 16 column file, 12 for the DOE factors, 4 for the PCR result types. Rows are the DOE test conditions followed by the respective 4 PCR results.
# sobol1.csv can also be a 12 column file with just the DOE factors if you dont instance.RUN()
```

Data Importing - Using a saved dataset:
```python
instance = DataBall()
instance.importFromDirectory()
~~~ input request: # select folder and FILE.csv through interactive menu.
~~~ input request: # Designate the number of input columns in the dataset (test factors).
~~~ input request: # Designate the number of output columns in teh dataset (results columns).
~~~ input request: # if ANALYTICS.mode() is on, it will ask if you want to override your temp data.

# creates a DataBall variable named: self.FILENAME_import<version> of type pd.DataFrame, and puts it into a cache.
# creates a self.factor_DF_FILENAME<version> file containing the information of the factors imported, and put its into a cache.
```

Data selection (IMPORTANT)
```python
instance.DOE_current_design() # shows what designs are currently available and activated to RUN(). CREATES self.DOE_active_pointer.
instance.DOE_current_design(change=int) # activates a different design to RUN(). UPDATES self.DOE_active_pointer FOR REST OF THE PROGRAM. 
instance.DOE_update(change=False) # returns min,max, current value for DOE factors. if change=True allows individual update of each. Empty string skips entry. 

# the program moves a pointer (self.DOE_active_pointer) to a local cache of all DOE designs, including the recently imported.
# the pointer will dictate what the Data calling, Data visualisation, Data Handling and Data Exporting functions will pull data from.
# the pointer will influence Temporary State Priming functions: it will dictate what dataset will be offered for temporary modifications. 
```

Data calling:
```python
instance.dataset() #returns imported DOE design or the full dataset (DOE_design + appended data) if RUN() has been executed.
instance.dataset_DOEmatrix() # returns just the DOE design.
isntance.dataset_results() # returns just the simulation data.
instance.dataset_reset() # resets heightened states and work done under them. Brings back the original DOE dataset.  
```

Data visualisation:
```python
instance.plot_byoutput() # Plots 4 scatter plots, 1 per output of the model. the x axis is the nr. of runs of your design.
instance.plot_byfactor(x_factors,regression_coefficients) # Plots 12 scatter plots, 1 per test factor. x axis is your factor, y axis is DNA yield only.
# regression coefficients can be generated through the ANOVA class and passed into the function to see the plotted linear regressions.
instance.plot_distributions(factor): 4 scatterplots: Normal Probability plot, Histogram, residuals by row, residuals by predicted. Data needs to be sorted for normal probability plot. Note that the last plot needs a predicted column, which can only be generated if you run linear regression.
instance.plot_withingroups(statistic=0): 1 plot: A comparison of a given statistic for an output result (DNA yield only) between each treatment levels of a test factor. In other words, a DOE plot. The DOE plot can be used to plot average results between treatments (statistic = 0) or standard deviations between treatments (statistic = 1).

PLOT STACKING: (instance.plot_byoutput(), .plot_distributions(factor), .plot_withingroups(statistic=0) 
These functions accept external data from other DataBall objects so that it can be overlayed into a single plot.
To stack plots, pass add=[DataBall2.dataset()]. For example DataBall1.plot_byoutput(add=[Databall2.dataset()]). Note that DataBall.ANALYTICS_mode() will mediate if the temp work data (if "ON") or cache data (if "OFF") will be used in EITHER of the Databall objects.

```

Data Handling:
```python
instance.data_topvalues(n=int,object=False) # Shows the top N values of the selected DOE matrix after RUN(). object=True allows assignment to a variable.
instance.data_sort(sort_ascending=True,object=False) # sorts data on the selected DOE matrix after RUN(). sort_ascending=False changes the sort, object=True allows assignment to a variable.
instance.data_average_ranking() # returns the average rank of each run across all 4 outputs as an appended column.
```

Data testing:
```python
ANOVA(DataBall_object,x:list,y:str): Generates an instance of ANOVA with the current version of the data from DataBall.dataset() and the factors meant for testing. x is a list with the input variable names. y is the output data, for example "DNA yield (ng/uL)".
  ANOVA_test() - useable when 1 x is passed (e.g. ["topprimer_vol"]. Generates a standard ANOVA results table.
  MANOVA_test() - useable when 2+ x are passed (e.g. ["topprimer_vol","dNTP_vol"]. Generates ANOVA results table.

  ANOVA.coefficient_matix() - generates a list of multiple linear regression factors (e.g. [b0,b1,b2,b3,b4...], where b0 is the intercept and b1... are the slopes of the respective factors).

  ~~ ANOVA() has all the underlying features necessary for the calculations of sum of squares, degrees of freedom, F_tests and R^2 calculations applicable to the regression models, error and total segments of the data, for both simple or muliple (via matrix multiplication). 

REGRESSION(DataBall,[x_list],y:str): Lets the user call for the estimation of linear regression factors, calculate AIC and BIC scores and apply stepwise minimisation of model paramerization (with AIC only)

MODEL_ANALYSIS(ANOVA_object): Lets the user create Y_predicted columns, residuals and standardized residuals for plotting. Must take an ANOVA object.
  MODEL_ANALYSIS.RUN(): adds Y_predicted and residuals for the list of factors passed to the ANOVA object.
  MODEL_ANALYSIS.std_resid(): adds standardized residuals.  
```

State priming:
```python
instance.ANALYTICS_mode()
# heightened state for working and visualising data. calling it will turn it ON or OFF. default = "OFF" 
# if OFF (default), all internal DataBall functions and other classes that request DataBall, will be editing the cache file directly.
# if ON, a temporary copy of the cache file is provided to the internal functions of DataBall and any other Class that requests it. All classes that request DataBall will edit this copy and a local copy of it simultaneously. 




