In order to build and run the code you have to use the BASH Unix shell.

Build: 
	To build the code you should navigate to hog2-PDB-refactor\build\gmake and run "make OPENGL=STUB" .

Run: 
	To run the tests that are presented in the paper you should navigate to hog2-PDB-refactor\bin\release and run "./bidirectional -aaai <filename>" , the filename parameter is the name you wish the resutls file will be named after. The default value is the current timestamp.
	This command will run small tests for the "Pancake Sorting" problem, the "Sliding Tile Puzzle" problem and grid problems.
	It is possible to change the parameters of the problems themselves in order to recreate the original tests of the paper.
	To modify the parameters navigate to hog2-PDB-refactor\apps\bidirectional\BidirTests.cpp file. 
	Pancake Sorting (lines 59-62): 
		There are 4 parameters - pancakes_num, randomPancakes, gaps, problem_num. set randomPancakes to False to get hard instances of pancake sorting, it is available only for pancakes_num = 16/20/24/28 and for maximum of 100 problems_num.
	Sliding Tile Puzzle (lines 100-102):
		There are 3 parameters - walkLength, randomSTP and problems_num. The walkLength used only when randomSTP=true. When randomSTP=false the problems used are known Korf difficult instances.
	Grid (lines 52+133):
		There are 2 parameters - mapName which defines on what grid the search is performed and problems_num.
	

