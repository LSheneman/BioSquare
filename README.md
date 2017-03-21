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
| 14       | Verbose Output                | T = ERH       |
| 15       | Hypothesis                    | F = BRH       |

Please note that Verbose Output produces the genomes of all species and population stats every 1023 generations.

The program can now be ran as:
	$ ./Biosquare {options}

Reproducting the Stats 

<<<<<<< HEAD
This script re-compiles data from the Biosquare output files. This R script can be run to reshape data and organize for analyses. All a statistical procedures, outputs and graphical results are generated from this R script.

First open the R file and change setwd to the path of your data files. The the R script can be run in either the GUI or in the command line using:

	$ Rscript PopPostRun.R

Note: Run these two scripts separately as they both utilize and reshape the data unique to the analysis.
=======
>>>>>>> cd31274af911724b0a4d8f8145e0e23f5897c893
