**How to use:**

From no prior data:
```python
instance = DataBall() # creates DataBall class object

instance.DOE_import(sobol) # imports a DOE design, in this case a sobol.

~~~ input request: # requests user input for design-specific variables. for sobol, it will ask for total nr. of runs.

instance.RUN() # runs the sobol design through the PCR simulator

instance.savetoDirectory() # generates a "year_month_date - Results Folder" in curdir and saves results as sobol1_data.csv
```

From saved data:
```python
instance = DataBall()

instance.importFromDirectory()

~~~ input request: # select folder and FILE.csv through interactive menu.

_creates a DataBall variable named: self.FILE_import<version>_

```

**Planned tasks:**
- append the DOE matrix to the output data for future ease of use when plotting data.
- Update and upgrade functions for DOE data visualisation to accept DataBall compatibility




**userguide:**

```python
DataBall()
# instantiating DataBall collects all factors and their data from the PCR simulator
# generates a DataBall.factor_DF DataFrame with the data.

DataBall.DOE_import(sobol)
 # calling .DOE_import(sobol) will call pyDOE.sobol_sequence(), passign as arguments those specified via of DataBall.function_mapping(sobol)
 # creates a DataBall.sobol<version> (working test ranges given Min Max values) and a DataBall.sobol<version>**a** variable (boilerplate matrix. factors are between 0 and 1)
 # caches DataBall1.sobol<version> and sobol<version>**a** into DataBall.DOE_cache, and moves DataBall.DOE_active_pointer to it.


DataBall.DOE_current_design(change="")
# shows you what's in the DOE cache and what DOE you have selected for testing.
# change = int changes which DOE is meant to be used for testing.

DataBall.RUN()
# Takes each row of the selected DOE, passes it into DataBall.factor_DF, formats it for the PCR simulation, runs the PCR simulation for all the rows in the selected DOE.
# **kwargs accepts {"hard_limit" = integer} to truncate the DOE design, for development reasons. 








