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


