# PyKinet
This script is used to fit the model by Ballauf et al. (https://doi.org/10.1021/jp101125j) to the kinetics data for 4-nitrophenol reduction with sodium borohydride catalyzed by metal nanoparticles 

# Guidelines for the input data for the script:
-- The directory with the script should contain two subdirectories named "raw" and "model"
-- The subdirectory "raw" contains files that a user wants to process
-- The subdirectory "model" contains files that the script produces while processing the corresponding files from "raw"
-- Names of files in "raw" should end with "_SS_raw.txt" or "_SS_raw.TXT" to be recognised by the script and get processed (SS stands for steady state).
   Also, the file name is expected to start with sX where X is one or several digits to denote the sample.
   Therefore, the file name should have the following format: "sXX(any characters)_SS_raw.txt"
-- The corresponding file with the results of processing is saved in "model" and it has the same name as the original data file but ends with "_model.txt"
-- The data files in "raw" should be formatted according to the template:
	- Block 1: parameter name and its value separated by a tab or a space	
		- The user doesn't have to set the parameters in the file but it's the best way to change them. If they are not set in the file, the default values from the script are used.
		- The parameters used:				
			Knp_init	-- the initial (before optimization) value of the Langmuir adsorption constant of 4-nitrophenol (a fitting parameter)		
			Kbh_init	-- the initial (before optimization) value of the Langmuir adsorption constant of sodium borohydride (a fitting parameter)
			k_init		-- the initial (before optimization) value of the rate constant k (a fitting parameter)
			n_init		-- the initial (before optimization) Langmuir−Freundlich exponent for 4-nitrophenol (a fitting parameter)
			m_init		-- the initial (before optimization) Langmuir−Freundlich exponent for sodium borohydride (a fitting parameter)
			Cbh			-- the concentration of sodium borohydride (a fixed parameter)
			S			-- the total surface area of nanoparticles (a fixed parameter)
		- According to the model, any of these parameters can be fixed or fitted 
		- Don't change the names of the parameters. Change only their values.
	- Separator line: 20 dash ('-') symbols . It is used to divide the Block 1 from the Block 2. This line is necessary for succeful parsing of the files
	- Block 2: the data starting with the line of headers
		- Starts with one line with the headers. It doesn't matter what the headers are
		- Should contain numbers in the rest of the lines. Numbers in different columns should be divided with tab or space characters.
		  This is applied automatically when you simply copy-paste data from Excel.

# Requirements
Executing the script requires installation of the following Python libraries: os, re, numpy, scipy, matplotlib, and lmfit.
The performance was tested only with Python 3.6.7, numpy 1.19.1, scipy 1.5.2, matplotlib 3.3.1, lmfit 1.0.1

# The script works in several steps:
1) Necessary libraries are imported, the default directories are set up
2) Default values for the parameters are set up
3) Looks for data files in the subdirectory /raw, the names of the files should end with _SS.txt to get recognized and processed
4) Reads the data files line by line and extracts the data according to the pattern of data files:
    a) Block 1: parameter name and its value separated by a tab or a space
    b) Separator line: 20 dash symbols
    c) Block 2: the data starting with the line of headers
5) Groups the data sets by samples (all replicats of the same sample belong to one group)
6) The system of ODEs and the objective function (the one that gets minimized by the solver) are described
7) For each data set fitting is set up and the optimal values of the parameters are found
8) Both input and fitted data are plotted
    A) For individual samples
    B) For grouped samples
9) Both input and fitted data are exported to a txt file
    A) For individual samples
    B) For grouped samples
10) All rate constants are exported to a txt file

# Contact
For any questions, please contact andrey.v.romanyuk@gmail.com
