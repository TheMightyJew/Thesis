#include "BidirTests.h"

#include "TestInfo.h"

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
#include <sstream>
#include <boost/format.hpp>
#include <ctime>
using namespace std;

static ofstream myfile;
const string RESULTS_DIR_PATH = "Test_Results/";
static void AAAI_Pancake(const string file, int problemsNum);
static void AAAI_STP(const string file, int problemsNum);
static void AAAI_Grid(const string file, int problemsNum, const char *mapName = "brc000d", double weight = 1.0);
static string getCurrentTime();

string getCurrentTime()
{
	time_t rawtime;
	struct tm *timeinfo;
	char buffer[80];

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer, 80, "%d-%m-%Y_%H-%M-%S", timeinfo);
	return string(buffer);
}

void AAAI_Test(string fileName)
{
	if (fileName.empty())
	{
		fileName = "results_" + getCurrentTime();
	}
	fileName += ".csv";

	const int problemsNum = 100;
	AAAI_Pancake(fileName, problemsNum);
	AAAI_STP(fileName, problemsNum);
	AAAI_Grid(fileName, problemsNum, "brc000d");

	exit(0);
}

void AAAI_Pancake(const string fileName, int problemsNum)
{
	bool randomPancake = true;
	vector<int> gaps = {0, 1, 2, 3};
	const int pancakesNum = 12;
	const string PROBLEM_NAME = "PancakeSorting";
	string filePath = RESULTS_DIR_PATH + PROBLEM_NAME + "/" + fileName;
	myfile.open(filePath);

	myfile << TestResult::csvSerializeHeaders() << std::endl;
	GenericTester<PancakePuzzleState<pancakesNum>, PancakePuzzleAction, PancakePuzzle<pancakesNum>, PancakeHasher<pancakesNum>> gt;
	for (int gap : gaps)
	{
		srandom(2017218);
		PancakePuzzleState<pancakesNum> original, goal;
		PancakePuzzle<pancakesNum> pancake(gap);
		pancake.SetUseRealValueEdges(false);

		std::stringstream testDescription;
		testDescription << boost::format("TestPancake:(Pancakes: %d, Gap: %d, Random: %d, Hard: %d)") % pancakesNum % gap % randomPancake % (!randomPancake);

		for (int count = 0; count < problemsNum; count++)
		{
			goal.Reset();
			original.Reset();
			if (randomPancake)
			{
				for (int x = 0; x < pancakesNum; x++)
					swap(original.puzzle[x], original.puzzle[x + random() % (pancakesNum - x)]);
			}
			else
			{
				GetPancakeInstance(original, count);
			}
			cout << "Running Pancake: Gap=" << gap << ", ProblemID=" << count + 1 << endl;

			std::stringstream startState;
			startState << original;
			std::stringstream goalState;
			goalState << goal;
			TestInfo testInfo = TestInfo(testDescription.str(), count, startState.str(), goalState.str());

			gt.genericTest(original, goal, pancake, myfile, testInfo);
		}
	}
	myfile.close();
}

void AAAI_STP(const string fileName, int problemsNum)
{
	srandom(2017218);
	int walkLength = 32;
	bool randomSTP = true;

	const string PROBLEM_NAME = "STP";
	string filePath = RESULTS_DIR_PATH + PROBLEM_NAME + "/" + fileName;
	myfile.open(filePath);
	myfile << TestResult::csvSerializeHeaders() << std::endl;
	GenericTester<MNPuzzleState<4, 4>, slideDir, MNPuzzle<4, 4>, STPHasher> gt;
	MNPuzzleState<4, 4> original, goal;
	MNPuzzle<4, 4> mnp;

	std::stringstream testDescription;
	testDescription << boost::format("TestSTP:(Random: %d, Hard: %d)") % randomSTP % (!randomSTP);

	for (int count = 0; count < problemsNum; count++)
	{
		goal.Reset();
		original.Reset();
		if (randomSTP)
		{
			original = STP::GetRandomInstance(walkLength);
		}
		else
		{
			original = STP::GetKorfInstance(count);
		}
		cout << "Running STP: ProblemID=" << count + 1 << endl;

		std::stringstream startState;
		startState << original;
		std::stringstream goalState;
		goalState << goal;
		TestInfo testInfo = TestInfo(testDescription.str(), count, startState.str(), goalState.str());

		gt.genericTest(original, goal, mnp, myfile, testInfo);
	}
	myfile.close();
}

void AAAI_Grid(const string fileName, int problemsNum, const char *mapName, double weight)
{
	const string PROBLEM_NAME = "Grid";
	string filePath = RESULTS_DIR_PATH + PROBLEM_NAME + "/" + fileName;
	myfile.open(filePath);
	myfile << TestResult::csvSerializeHeaders() << std::endl;
	GenericTester<CanonicalGrid::xyLoc, CanonicalGrid::tDirection, CanonicalGrid::CanonicalGrid, GridHasher> gt;
	CanonicalGrid::xyLoc original, goal;
	CanonicalGrid::CanonicalGrid *cg_tmp = 0;

	std::stringstream testDescription;
	testDescription << boost::format("TestGrid:( Map_Name: %s)") % mapName;

	srandom(2017218);
	MapEnvironment *me = 0;

	const char *scen_start_path = "../../scenarios/dao/";
	const char *scen_end_path = ".map.scen";
	string total_scen(string(scen_start_path) + mapName + scen_end_path);
	const char *scen_path = total_scen.c_str();

	const char *map_start_path = "../../maps/dao/";
	const char *map_end_path = ".map";
	string total_map(string(map_start_path) + mapName + map_end_path);
	const char *map_path = total_map.c_str();

	ScenarioLoader s(scen_path);
	Map *m = new Map(map_path);
	me = new MapEnvironment(m);
	me->SetDiagonalCost(1.5);
	cg_tmp = new CanonicalGrid::CanonicalGrid(m);
	CanonicalGrid::CanonicalGrid cg = *cg_tmp;
	cg.SetDiagonalCost(1.5);
	int experiments_num = s.GetNumExperiments();
	int counter = 0;
	for (int count = 0; count < experiments_num && counter < problemsNum; count++)
	{
		Experiment e = s.GetNthExperiment(count);
		if (e.GetDistance() <= 0)
			continue;

		original.x = e.GetStartX();
		original.y = e.GetStartY();
		goal.x = e.GetGoalX();
		goal.y = e.GetGoalY();

		cout << "Running Grid: ProblemID=" << counter + 1 << endl;

		std::stringstream startState;
		startState << original;
		std::stringstream goalState;
		goalState << goal;
		TestInfo testInfo = TestInfo(testDescription.str(), count, startState.str(), goalState.str());

		gt.genericTest(original, goal, cg, myfile, testInfo);
		counter++;
	}
	myfile.close();
}