Data replication:

Set up BioSquare:

Download the source code by executing the following for commands:
	
	$ git clone https://github.com/LSheneman/BioSquare.git
	$ cd BioSquare
	$ git init

Next you will need to compile the code.
	
	$ g++ -std=c++11 *.cpp -O3 -o Biosquare


Defining Values and Running the Code:

The following is an order of the commandline arugements for BioSquare. They must be put in order. If there is no argument passed in the default value is used:


| Argument | Option                        | Default Value |
|----------|-------------------------------|---------------|
| 1        | Population file name          | NULL          |
| 2        | Mutation Rate                 | 0.01          |
| 3        | Replacement Rate              | 0.01          |
| 4        | Herbivore to Plant Payoff     | -1.0          |
| 5        | Microbes to Plant Payoff      | 1.0           |
| 6        | Plant to Herbivore Payoff     | 1.0           |
| 7        | Plant to Microbes Payoff      | 1.0           |
| 8        | Plant X Starting Location     | 0             |
| 9        | Plant Y Starting Location     | 0             |
| 10       | Herbivore X Starting Location | 0             |
| 11       | Herbivore Y Starting Location | 0             |
| 12       | Microbes X Starting Location  | 0             |
| 13       | Microbes Y Starting Location  | 0             |
| 14       | Verbose Output                | False         |
| 15       | Biotic Resistance Hypothesis  | False         |

The program can now be ran as:
	$ ./Biosquare {options}


Data Preparation for Analyses with R Scripts:


Optimal Search Statistics 

Once the local directory of raw Avida data is stored as text files this R script can be run to reshape data and organize for analyses. 
All a statistical procedures, outputs and graphical results pertaining to figures 1 ("Average Final Gestation Time by World Size"), 2 ("Average Genome Length by World Size"), and 4 ("Average Gestion Time by Genome Length") are generated from this R script.

In the context of this tutorial, the "geno_XX.txt" files can be found in the docs/data directory. To run the scripts from the command line we must first place them in the same folder as the data.
	
	$ cd docs/data
	$ cp ../../scripts/* .

You can run the script to obtain optimal search statistics:
	$ Rscript search_stats.R

Instruction Propertions Stats 

This script re-compiles data from the Avida output files and the using the genome text files already saved in local directories (see above). This R script can be run to reshape data and organize for analyses. All a statistical procedures, outputs and graphical results pertaining to figure 3 ("Propertions of Non-Movement Instructions") is generated from this R script.

To run it, ensure that you are in the docs/data folder.
	$ pwd 

will output your current directory. Make sure this path ends with "/docs/data".

Now you can produce a visual of the instruction propertion stats by typing:
	$ Rscript proportion_stat.R

Note: Run these two scripts separately as they both utilize and reshape the data unique to the analysis.