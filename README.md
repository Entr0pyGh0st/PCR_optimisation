Experimental designs are methods of querying the behaviour of complex in silico/in vivo systems in an extensive and data-efficient manner.

They address a technical execution problem that arises when a scientist/engineer seeks to investigate and model a black box system with large sets of input parameters empirically. Generally speaking, both the inputs and output are usually continous and the resulting approximation model is usually a simple combination of multiple linear or quadratic models. 

Whilst individual input parameter testing is feasible, accepteable and of straight forward analysis and experimental control, it's cumbersome nature of iterating through each individual factor is both resource ineficient and "low resolution" (i.e. doesn't explain how inputs can interact with eachother to impact the output, only how each input affects the output exclusively).   

A subset of experimental designs addresses this by providing bespoke testing schematics that leverage controlled input aliasing (i.e. having two factors changing simultaneously in a test run) to 1) to test more inputs per test run and 2) investigate for input interactions (n-factor interactions), thus allowing for a more complete description of black box behaviour compared straightforward indididual factor testing. 

Another subset of experimental designs use quasirandom sequences of numbers (low-discrepancy sequences) to generate a testing schematic that covers the experimental ranges of the input factors of a black box evenly. These are technically superior to the latter subset if not constrained by resources or time.

Ben Shirt-Ediss made https://virtual-pcr.ico2s.org/pcr/, an in silico model of a PCR reaction aimed at amplifying a 1kb DNA sequence by changing the parameterization of 12 different inputs. The source code can be found at https://bitbucket.org/ben_s_e/virtual-pcr-notebook/src/main/ .Understandeably, the fact that he coded a model means that one could derive the maximum output metrics (yield, product purity) numerically. However, he and I (as found in this repo), intended that the solution be found via experimental designs, for the sakes of learning about them.

This repo seeks to interface his model with pyDOE, a python module for experimental designs, to provide an interactive approach to experimental design testing without going through the web server. 

------------------------------------------------------------------------------------------

## Current highscore: 1.066 mg/mL, 99.6% pure, 330 minute run - 1024 run SOBOL sequence

------------------------------------------------------------------------------------------

**How to use:**

- Create a DataBall instance.

- Call DataBall.DOE_import(sobol) -> imports a sobol design

- input: number of runs for sobol design

- Call DataBall.RUN()


**userguide:**

instantiating DataBall collects all factors and their data from the PCR simulator, whilst also generating a DataBall.factor_DF DataFrame for organisation purposes.

DataBall.DOE_import(sobol) runs the sobol design function with the required parameters within DataBall. DataBall.Function_Mapping is called to map the function to it's function call parameters.

DataBall.DOE_import(sobol) creates a DataBall.sobol1 (test levels updated) and a DataBall.sobol1a variable (boilerplate matrix)

DataBall.DOE_import(sobol) also caches DataBall.sobol1.

Recursively calling the prior function creates and caches DataBall.sobol2, sobol3, sobol4 etc... You can customise the sobol design with each call (currently it's just run nr.)

DataBall.DOE_current_design() shows you what's in the DOE cache and what DOE you have selected for testing. 

DataBall.DOE_current_design(change=int) changes the active DOE design to one from the cache.

DataBall.RUN() Takes each row of the selected DOE, passes it into DataBall.factor_DF, formats it for the PCR simulation, runs the PCR simulation for all the rows in the selected DOE. **kwargs accepts {"hard_limit" = integer} to truncate the DOE design, for development reasons. 


**Planned tasks:**
- Migrate functions for DOE run data collection into DataBall
- Update and upgrade functions for DOE data visualisation to accept DataBall compatibility
- Update and upgrade functions for locally saving data to accept DataBall compatibility






