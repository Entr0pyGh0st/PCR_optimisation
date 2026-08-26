**How to use:**

Data Generation - Generating a fresh dataset:
```python
instance = DataBall() # creates DataBall class object
instance.DOE_import(sobol) # imports a DOE design, in this case from the function named sobol. Creates an internal attribute named self.<design><version>.
~~~ input request: # requests user input for design-specific variables. for sobol, it will ask for total nr. of runs.
instance.RUN() # runs the newly created attribute ,e.g. self.sobol1,  through the PCR simulator

```

Data Exporting - Saving a dataset:
```python
instance.savetoDirectory() # generates a "year_month_date - Results Folder" in curdir and saves results as sobol1.csv

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

# the program moves a pointer (self.DOE_active_pointer) to a local cache of all DOE designs, including the recently imported.
# the pointer will dictate what the Data calling, Data visualisation and Data Handling funcitons will pull data from.
# the pointer will influence Temporary State Priming functions: it will dictate what dataset will be offered for temporary modifications. 
```

Data calling:
```python
instance.dataset() #returns imported DOE design or the full dataset (DOE_design + appended data) if RUN() has been executed.
instance.dataset_DOEmatrix() # returns just the DOE design.
isntance.dataset_results() # returns just the simulation data.
instance.dataset_reset() # resets heightened states and work done unter them. resets to the DOE design selected with the Data Selection functions.  
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
```python:
instance.ANALYTICS_mode()
# heightened state for working and visualising data. calling it will turn it ON or OFF. default = "OFF" 
# if OFF (default), ALL Data Handling functions WILL RETURN AN OBJECT BY DEFAULT, overriding object=True assignments. 
# if ON, ALL Data Calling, Data Handling, Data visualisation and Data Saving functions will now be applicable to an inplace copy of active DOE file. This means that all Data Handling functions now behave as object=False by default (don't allow for variable assignment).
# However, you can still pass object=True when ANALYTICS_mode = ON to assign a variable to the result of a Data Handling functio.
```

**Planned tasks:**
- Add functionality to the functions in Data Visualisation to accept an external variable instead of built-in data calling methods.
- ANALYTICS_mode() clashes with Data Importing and Data Selection because of forced self.DOE_active reassignment, thus deleting the inplace copy of the active DOE file.






