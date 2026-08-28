

Experimental designs are methods of querying the behaviour of complex in silico/in vivo systems in an extensive and data-efficient manner.

They address a technical execution problem that arises when a scientist/engineer seeks to investigate and model a black box system with large sets of input parameters empirically. Generally speaking, both the inputs and output are usually continous and the resulting approximation model is usually a simple combination of multiple linear or quadratic models. 

Whilst individual input parameter testing is feasible, accepteable and of straight forward analysis and experimental control, it's cumbersome nature of iterating through each individual factor is both resource ineficient and "low resolution" (i.e. doesn't explain how inputs can interact with eachother to impact the output, only how each input affects the output exclusively).   

A subset of experimental designs addresses this by providing bespoke testing schematics that leverage controlled input aliasing (i.e. having two factors changing simultaneously in a test run) to 1) to test more inputs per test run and 2) investigate for input interactions (n-factor interactions), thus allowing for a more complete description of black box behaviour compared straightforward indididual factor testing. 

Another subset of experimental designs use quasirandom sequences of numbers (low-discrepancy sequences) to generate a testing schematic that covers the experimental ranges of the input factors of a black box evenly. These are technically superior to the latter subset if not constrained by resources or time.

Ben Shirt-Ediss made https://virtual-pcr.ico2s.org/pcr/, an in silico model of a PCR reaction aimed at amplifying a 1kb DNA sequence by changing the parameterization of 12 different inputs. The source code can be found at https://bitbucket.org/ben_s_e/virtual-pcr-notebook/src/main/ .Understandeably, the fact that he coded a model means that one could derive the maximum output metrics (yield, product purity) numerically. However, he and I (as found in this repo), intended that the solution be found via experimental designs, for the sakes of learning about them.

This repo seeks to interface his model with pyDOE, a python module for experimental designs, to provide an interactive approach to experimental design testing without going through the web server. 

------------------------------------------------------------------------------------------

## Current highscore: 1.483 mg/mL, 97.6% pure, 8801-fold amplification, 22060 second run

------------------------------------------------------------------------------------------

Bugs:
-instance.importFromDirectory() doesn't revert back to original directory, thus recursively saving files.
-instance.importFromDirectory() forces deletion of temp data generated under ANALYTICS_mode()
-instance.importFromDirectory() breaks when importing data with Average Rankings column appended. 
-instance.DOE_current_design() will override temp file work enabled by instance.ANALYTICS_mode()

future features:
-let Data Visualisation functions use external data.
-statistical testing for interactions
-stepwise minimisation of BIC and AIC(c) of models.
-more designs: Fractional, PB, BB
-more plots: DOE mean, stdev, scatter and interaction plots ,residual by row, residual by predicted, Q-Q,

**How to use:**

Data Generation - Generating a fresh dataset:
```python
instance = DataBall() # creates DataBall class object
instance.DOE_import(sobol) # imports a DOE design, in this case from the function named sobol. Creates an internal attribute named self.<design><version>.
~~~ input request: # requests user input for design-specific variables. for sobol, it will ask for total nr. of runs.
instance.RUN() # runs the newly created attribute ,e.g. self.sobol1,  through the PCR simulator

DOE_update(True) followed by DOE_import(sobol) # 1st function call updates the min,max values. 2nd function call generates an updated sobol design.

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

# creates a DataBall variable named: self.FILE_import<version> of type pd.DataFrame
```

Data selection (IMPORTANT)
```python
instance.DOE_current_design() # shows what designs are currently available and activated to RUN()
instance.DOE_current_design(change=int) # activates a different design to RUN()
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
instance.plot() # 4 scatterplots, 1 per PCR result type. (y = PCR result // x = test run nr.)
instance.plot_byfactor() # 12 scatterplots, 1 per factor. (y = PCR result // x = factor values) ## ONLY PLOTS DNA YIELD CURRENTLY.
```

Data Handling:
```python
instance.data_topvalues(n=int,object=False) # Shows the top N values of the selected DOE matrix after RUN(). object=True allows assignment to a variable.
instance.data_sort(sort_ascending=True,object=False) # sorts data on the selected DOE matrix after RUN(). sort_ascending=False changes the sort, object=True allows assignment to a variable. 
```

State priming:
```python
instance.ANALYTICS_mode()
# heightened state for working and visualising data. calling it will turn it ON or OFF. default = "OFF" 
# if OFF (default), ALL Data Handling functions WILL RETURN AN OBJECT BY DEFAULT. the underlying dataset will not be modified.
# if ON, ALL Data Calling, Data Handling, Data visualisation and Data Saving functions will now be applicable to an inplace copy of the underlying dataset. This means that all Data Handling functions now behave as object=False.
# However, you can still pass object=True when ANALYTICS_mode = ON to assign a variable to any modification done to the working copy of the dataset.
```




