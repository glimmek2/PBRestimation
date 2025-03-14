# PBRestimation
•	calling_program_for_all_scenarios.R is a wrapper for all simulations of the paper in one go.
• makeplot_Fig3.R generates Figure 3 of the paper. 

•	calling_program_for_simulation_runs.R is the «master code” that distributes the simulations of a single scenario onto many cores.
•	oneSimRun.R is the “submaster code” which calls the different functions for first simulating, then analyzing the data.
•	simulation_functions_myseed.R contains all the functions needed to simulate data.
•	analysis_functions_sep_robustvar_fix.R contains all the functions needed to analyze data.

• makeplot.R is a sample call generating a figure of PBR curves and difference curve.
• eventsdata.csv is the dataset for figure 3.
