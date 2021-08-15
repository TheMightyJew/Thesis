#include "BidirTests.h"

#include "TestInfo.h"

#include "MNPuzzle.h"
#include "STPInstances.h"

#include "PancakePuzzle.h"
#include "PancakeInstances.h"

#include "CanonicalGrid.h"
#include "ScenarioLoader.h"
#include "Map2DEnvironment.h"
#include "MapGenerators.h"
#include "MapOverlay.h"
#include "json.h"

#include "GenericTester.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <boost/format.hpp>
#include <sys/stat.h>
#include <sys/types.h>
#include <ctime>

using namespace std;
using json = nlohmann::json;

static ofstream myfile;
const string OUTPUT_DIR_PATH = "BidirTests/";
const string RESULTS_DIR_PATH = OUTPUT_DIR_PATH + "Results/";
const string CONFIG_FILENAME = "config.json";
const string CONFIG_PATH = OUTPUT_DIR_PATH + CONFIG_FILENAME;
const string DEFAULT_CONFIG_PATH = OUTPUT_DIR_PATH + "default_config.json";
static void AAAI_Pancake(const string test_directory_path, json configurations);
static void AAAI_STP(const string test_directory_path, json configurations);
static void AAAI_Grid(const string test_directory_path, json configurations);
static string getCurrentTime();

string getCurrentTime()
{
	const int bufferSize = 80;
	const char* timeFormat = "%d-%m-%Y_%H-%M-%S";
	time_t rawtime;
	struct tm *timeinfo;
	char buffer[bufferSize];

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer, bufferSize, timeFormat, timeinfo);
	return string(buffer);
}

void AAAI_Test(string test_name)
{
	if (test_name.empty())
	{
		test_name = "test_" + getCurrentTime();
	}
	string test_directory_path = RESULTS_DIR_PATH + test_name;
	mkdir(test_directory_path.c_str(), 0777);

	json configurations;
	ifstream config_file(CONFIG_PATH);
	ofstream copied_config_file(test_directory_path + "/" + CONFIG_FILENAME);
	if (!config_file.is_open())
	{
		config_file = ifstream(DEFAULT_CONFIG_PATH);
	}
	config_file >> configurations;
	string line;
	config_file.clear();
	config_file.seekg(0);
	while (getline(config_file, line))
	{
		copied_config_file << line << "\n";
	}
	config_file.close();
	copied_config_file.close();

	if (configurations["problems"]["pancake"]["shouldRun"].get<bool>())
	{
		AAAI_Pancake(test_directory_path, configurations);
	}
	if (configurations["problems"]["stp"]["shouldRun"].get<bool>())
	{
		AAAI_STP(test_directory_path, configurations);
	}
	if (configurations["problems"]["grid"]["shouldRun"].get<bool>())
	{
		AAAI_Grid(test_directory_path, configurations);
	}

	exit(0);
}

void AAAI_Pancake(const string test_directory_path, json configurations)
{
	bool randomPancake = !configurations["problems"]["pancake"]["hard"].get<bool>();
	vector<int> gaps = configurations["problems"]["pancake"]["gaps"];
	const int pancakesNum = 10;
	//int pancakesNum =  configurations["problems"]["pancake"]["size"].get<int>();
	const string PROBLEM_NAME = "PancakeSorting";
	string filePath = test_directory_path + "/" + PROBLEM_NAME + ".csv";
	myfile.open(filePath);

	myfile << TestResult::csvSerializeHeaders() << std::endl;
	GenericTester<PancakePuzzleState<pancakesNum>, PancakePuzzleAction, PancakePuzzle<pancakesNum>> gt;
	for (int gap : gaps)
	{
		srandom(2017218);
		PancakePuzzleState<pancakesNum> original, goal;
		PancakePuzzle<pancakesNum> pancake(gap);
		pancake.SetUseRealValueEdges(false);

		std::stringstream testDescription;
		testDescription << boost::format("TestPancake:(Pancakes: %d, Gap: %d, Hard: %d)") % pancakesNum % gap % (!randomPancake);

		for (int count = 0; count < configurations["problems"]["quantity"].get<int>(); count++)
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
			if (count == 7 || count == 39)
			{
				cout << "Skipped Pancake: Gap=" << gap << ", ProblemID=" << count + 1 << endl;
				continue;
			}
			cout << "Running Pancake: Gap=" << gap << ", ProblemID=" << count + 1 << endl;

			std::stringstream startState;
			startState << original;
			std::stringstream goalState;
			goalState << goal;
			TestInfo testInfo = TestInfo(testDescription.str(), count, startState.str(), goalState.str());

			gt.genericTest(original, goal, pancake, myfile, testInfo, 1, configurations);
		}
	}
	myfile.close();
}

void AAAI_STP(const string test_directory_path, json configurations)
{
	srandom(2017218);
	const int WALK_LENGTH = 32;
	const int SIZE = 4;
	bool randomSTP = !configurations["problems"]["stp"]["hard"].get<bool>();

	const string PROBLEM_NAME = "STP";
	string filePath = test_directory_path + "/" + PROBLEM_NAME + ".csv";
	myfile.open(filePath);
	myfile << TestResult::csvSerializeHeaders() << std::endl;
	GenericTester<MNPuzzleState<SIZE, SIZE>, slideDir, MNPuzzle<SIZE, SIZE>> gt;
	MNPuzzleState<SIZE, SIZE> original, goal;
	MNPuzzle<SIZE, SIZE> mnp;

	std::stringstream testDescription;
	testDescription << boost::format("TestSTP:(Hard: %d)") % (!randomSTP);

	for (int count = 0; count < configurations["problems"]["quantity"].get<int>(); count++)
	{
		goal.Reset();
		original.Reset();
		if (randomSTP)
		{
			original = STP::GetRandomInstance(WALK_LENGTH);
		}
		else
		{
			original = STP::GetKorfInstance(count);
		}

		if (original == goal)
		{
			continue;
		}
		cout << "Running STP: ProblemID=" << count + 1 << endl;
		std::stringstream startState;
		startState << original;
		std::stringstream goalState;
		goalState << goal;
		TestInfo testInfo = TestInfo(testDescription.str(), count, startState.str(), goalState.str());

		gt.genericTest(original, goal, mnp, myfile, testInfo, 1, configurations);
	}
	myfile.close();
}

void AAAI_Grid(const string test_directory_path, json configurations)
{
	const string PROBLEM_NAME = "Grid";
	string filePath = test_directory_path + "/" + PROBLEM_NAME + ".csv";
	myfile.open(filePath);
	myfile << TestResult::csvSerializeHeaders() << std::endl;
	GenericTester<CanonicalGrid::xyLoc, CanonicalGrid::tDirection, CanonicalGrid::CanonicalGrid> gt;
	CanonicalGrid::xyLoc original, goal;
	CanonicalGrid::CanonicalGrid *cg_tmp = 0;
	string mapName = configurations["problems"]["grid"]["map"].get<string>();
	cout << mapName << endl;
	std::stringstream testDescription;
	testDescription << boost::format("TestGrid:( Map_Name: %s)") % mapName;
	cout << mapName << endl;
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
	for (int count = 0; count < experiments_num && counter < configurations["problems"]["quantity"].get<int>(); count++)
	{
		Experiment e = s.GetNthExperiment(count);
		if (e.GetDistance() <= 0)
		{
			continue;
		}
		original.x = e.GetStartX();
		original.y = e.GetStartY();
		goal.x = e.GetGoalX();
		goal.y = e.GetGoalY();

		cout << "Running Grid: ProblemID=" << counter + 1 << endl;

		std::stringstream startState;
		startState << original;
		std::stringstream goalState;
		goalState << goal;
		TestInfo testInfo = TestInfo(testDescription.str(), counter, startState.str(), goalState.str());

		gt.genericTest(original, goal, cg, myfile, testInfo, 0.5, configurations);
		counter++;
	}
	myfile.close();
}