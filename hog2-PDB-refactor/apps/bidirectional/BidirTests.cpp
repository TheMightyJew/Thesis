#include "BidirTests.h"

#include "BidirSTP.h"
#include "MNPuzzle.h"
#include "STPInstances.h"
#include "STPHasher.h"

#include "PancakePuzzle.h"
#include "PancakeInstances.h"
#include "PancakeHasher.h"

#include "CanonicalGrid.h"
#include "ScenarioLoader.h"
#include "Map2DEnvironment.h"
#include "MapGenerators.h"
#include "MapOverlay.h"
#include "GridHasher.h"

#include "GenericTester.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <boost/format.hpp>
#include <ctime>
using namespace std;



static ofstream myfile;
static string filename;
static void AAAI_Pancake(string file);
static void AAAI_STP(string file);
static void AAAI_Grid(string file, const char *mapName = "brc000d", double weight = 1.0);
static string datetime();

string datetime()
{
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time (&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer,80,"%d-%m-%Y_%H-%M-%S",timeinfo);
    return string(buffer);
}

void AAAI_Test(string file){
	
	AAAI_Pancake(file);
	AAAI_STP(file);
	AAAI_Grid(file);
	
	exit(0);
}


void AAAI_Pancake(string file){
	const int pancakes_num = 8;
	bool randomPancake = true;
	vector<int> gaps = {0, 1, 2};
	int problems_num = 10;
	
	if(file == "")
		file = "results_" + datetime();
	filename = "Test_Results/PancakeSorting/" + file + ".txt";
	cout << "running..." << endl;
	myfile.open(filename);
	
	GenericTester<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, PancakePuzzle<pancakes_num>, PancakeHasher<pancakes_num>> gt;
	for(int gap : gaps){
		srandom(2017218);
		PancakePuzzleState<pancakes_num> original, goal;
		PancakePuzzle<pancakes_num> pancake(gap);
		pancake.SetUseRealValueEdges(false);
		
		myfile << boost::format("TestPancake:(Pancakes: %d, Gap: %d, Random: %d, Hard: %d)\n") % pancakes_num % gap % randomPancake % (!randomPancake);
		for (int count = 0; count < problems_num; count++){
			goal.Reset();
			original.Reset();
			if(randomPancake){
				for (int x = 0; x < pancakes_num; x++)
					swap(original.puzzle[x], original.puzzle[x+random()%(pancakes_num-x)]);
			}
			else{
				GetPancakeInstance(original, count);
			}
			cout << "Running Pancake: Gap=" << gap << ", ProblemID=" << count+1 << endl;
			myfile << boost::format("\tProblem %d of %d\n") % (count+1) % problems_num;
			gt.genericTest(original, goal, pancake,  myfile);
		}
	}
	myfile << "completed!" << endl;
	myfile.close();
	cout << "completed!" << endl;
}

void AAAI_STP(string file){
	srandom(2017218);
	int walkLength = 32;
	bool randomSTP = true;
	int problems_num = 10;
	
	if(file == "")
		file = "results_" + datetime();
	filename = "Test_Results/STP/" + file + ".txt";
	cout << "running..." << endl;
	myfile.open (filename);
	
	GenericTester<MNPuzzleState<4, 4>, slideDir, MNPuzzle<4, 4>, STPHasher> gt;
	MNPuzzleState<4, 4> original, goal;
	MNPuzzle<4, 4> mnp;
	myfile << boost::format("TestSTP:(Random: %d, Hard: %d)\n") % randomSTP % (!randomSTP);
	for (int count = 0; count < problems_num; count++){
		goal.Reset();
		original.Reset();
		if(randomSTP){
			original = STP::GetRandomInstance(walkLength);
		}
		else{
			original = STP::GetKorfInstance(count);
		}
		cout << "Running STP: ProblemID=" << count+1 << endl;
		myfile << boost::format("\tProblem %d of %d\n") % (count+1) % problems_num;
		gt.genericTest(original, goal, mnp, myfile);
	}
	myfile << "completed!" << endl;
	myfile.close();
	cout << "completed!" << endl;
}

void AAAI_Grid(string file, const char *mapName, double weight){
	if(file == "")
		file = "results_" + datetime();
	filename = "Test_Results/Grid/" + file + ".txt";
	cout << "running..." << endl;
	myfile.open (filename);
	GenericTester<CanonicalGrid::xyLoc, CanonicalGrid::tDirection, CanonicalGrid::CanonicalGrid, GridHasher> gt;
	CanonicalGrid::xyLoc original, goal;
	CanonicalGrid::CanonicalGrid *cg_tmp = 0;
	myfile << boost::format("TestGrid:( Map_Name: %s)\n") % mapName;
	srandom(2017218);
	MapEnvironment *me = 0;
	
	const char *scen_start_path = "../../scenarios/dao/";
	const char *scen_end_path = ".map.scen";
	string total_scen( string(scen_start_path) + mapName + scen_end_path);
	const char *scen_path = total_scen.c_str();

	const char *map_start_path = "../../maps/dao/";
	const char *map_end_path = ".map";
	string total_map( string(map_start_path) + mapName + map_end_path);
	const char *map_path = total_map.c_str();
	
	ScenarioLoader s(scen_path);
	Map *m = new Map(map_path);
	me = new MapEnvironment(m);
	me->SetDiagonalCost(1.5);
	cg_tmp = new CanonicalGrid::CanonicalGrid(m);
	CanonicalGrid::CanonicalGrid cg = *cg_tmp;
	cg.SetDiagonalCost(1.5);
	int problems_num = 10;
	int experiments_num = s.GetNumExperiments();
	int counter = 0;
	for (int count = 0; count < experiments_num && counter < problems_num; count++)
	{
		Experiment e = s.GetNthExperiment(count);
		if (e.GetDistance() <= 0)
			continue;
		
		original.x = e.GetStartX();
		original.y = e.GetStartY();
		goal.x = e.GetGoalX();
		goal.y = e.GetGoalY();
		
		cout << "Running Grid: ProblemID=" << counter+1 << endl;
		myfile << boost::format("\tProblem %d of %d\n") % (counter+1) % min(problems_num, experiments_num);
		gt.genericTest(original, goal, cg, myfile);
		counter++;
	}
	myfile << "completed!" << endl;
	myfile.close();
	cout << "completed!" << endl;
}