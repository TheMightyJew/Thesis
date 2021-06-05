/*
 *  GenericTester.hpp
 *
 *  Created by Steven Danishevski on SEP 20.
 *
 */

#ifndef GenericTester_H
#define GenericTester_H

#include <climits>
#include <algorithm>
#include "TemplateAStar.h"
#include "AStarOpenClosed.h"
#include "GenericTester.h"
#include "IDAStar.h"
#include "MM.h"
#include "MBBDS.h"
#include "BFBDS.h"
#include "IDBiHS.h"
#include "IDTHSwTrans.h"
#include "MbbdsBloomFilter.h"
#include "TestInfo.h"
#include "TestResult.h"
#include <fstream>
#include <iostream>
#include <boost/format.hpp>
#include <math.h>
#include <ctime>
using namespace std;

template <class state, class action, class environment, class hasher>
class GenericTester
{
private:
	const unsigned long STATES_QUANTITY_BOUND_DEFAULT = 1000000;
	const double MINIMUM_STATES_QUANTITY_BOUND = 2;
	const int SECONDS_LIMIT = 60 * 30;
	const vector<double> RMA_STATES_PERCENTAGES{0.5, 0.1, 0.01};

	bool AstarRun = true;
	bool RevAstarRun = true;
	bool IDAstarRun = true;
	bool AstarPIDAstarRun = true;
	bool AstarPIDAstarReverseRun = true;
	bool AstarPIDAstarReverseMinHRun = true;
	bool IDTHSpTrans = true;
	bool BAI = true;
	bool Max_BAI = true;
	bool MMRun = true;

	bool IDBiHSRun = true;
	bool F2Fheuristics = true;

	bool ASTARpIDBiHS = true;

	bool BFBDSRUN = true;
	bool revAlgo = false;
	bool threePhase = true;

	bool detectDuplicate = true;
	bool isConsistent = true;
	bool isUpdateByWorkload = true;
	bool printAbstractAlgos = false;

public:
	GenericTester() {}
	virtual ~GenericTester() {}
	void genericTest(state original, state goal, environment env, ofstream &myfile, TestInfo testInfo);
	vector<TestResult> testUMA(state original, state goal, environment env, TestInfo testInfo);
	vector<TestResult> testLMA(state original, state goal, environment env, TestInfo testInfo);
	vector<TestResult> testFMA(state original, state goal, environment env, TestInfo testInfo, unsigned long statesQuantityBound, vector<double> quantityPercentages);
};

template <class state, class action, class environment, class hasher>
void GenericTester<state, action, environment, hasher>::genericTest(state original, state goal, environment env, ofstream &myfile, TestInfo testInfo)
{
	vector<TestResult> lmaResults = testLMA(original, goal, env, testInfo);

	vector<TestResult> umaResults = testUMA(original, goal, env, testInfo);
	// Determine states quantity bound
	unsigned long statesQuantityBound = ULONG_MAX;
	for (TestResult testResult : umaResults)
	{
		if (0 < testResult.m_maxStatesInMemory && testResult.m_maxStatesInMemory < statesQuantityBound)
		{
			statesQuantityBound = testResult.m_maxStatesInMemory;
		}
	}
	if (statesQuantityBound == ULONG_MAX)
	{
		statesQuantityBound = STATES_QUANTITY_BOUND_DEFAULT;
	}
	vector<TestResult> fmaResults = testFMA(original, goal, env, testInfo, statesQuantityBound, RMA_STATES_PERCENTAGES);

	vector<vector<TestResult>> testResultsVectorts{lmaResults, umaResults, fmaResults};
	for (vector<TestResult> testResultsVector : testResultsVectorts)
	{
		for (TestResult testResult : testResultsVector)
		{
			myfile << testResult.csvSerialize() << std::endl;
		}
	}
}

template <class state, class action, class environment, class hasher>
vector<TestResult> GenericTester<state, action, environment, hasher>::testUMA(state original, state goal, environment env, TestInfo testInfo)
{
	double initialHeuristic = env.HCost(original, goal);
	Timer timer;
	vector<state> astarPath;

	vector<TestResult> testResults;
	// A*
	vector<AStarOpenClosedDataWithF<state>> astarOpenList;
	if (AstarRun)
	{
		TestResult astarTestResult = TestResult(testInfo);
		astarTestResult.m_initialHeuristic = initialHeuristic;
		astarTestResult.m_algorithmInfo = "A*";

		TemplateAStar<state, action, environment> astar;
		state start = original;
		timer.StartTimer();
		bool solved = astar.GetPathTime(&env, start, goal, astarPath, SECONDS_LIMIT);
		timer.EndTimer();
		if (solved)
		{
			astarTestResult.m_solutionCost = env.GetPathLength(astarPath);
			astarTestResult.m_statesExpanded = astar.GetNodesExpanded();
			astarTestResult.m_neccessaryStatesExpanded = astar.GetNecessaryExpansions();
			astarTestResult.m_maxStatesInMemory = astar.getMemoryStatesUse();
		}
		astarTestResult.m_timeElapsed = timer.GetElapsedTime();
		testResults.push_back(astarTestResult);
		if (printAbstractAlgos)
		{
			TestResult iastarTestResult = TestResult(testInfo);
			iastarTestResult.m_initialHeuristic = astarTestResult.m_initialHeuristic;
			iastarTestResult.m_algorithmInfo = "I-A*";
			iastarTestResult.m_solutionCost = astarTestResult.m_solutionCost;
			iastarTestResult.m_statesExpanded = astar.getIAstarExpansions();
			iastarTestResult.m_timeElapsed = astarTestResult.m_timeElapsed;
			testResults.push_back(iastarTestResult);
		}
	}
	if (RevAstarRun)
	{
		TestResult revAstarTestResult = TestResult(testInfo);
		revAstarTestResult.m_initialHeuristic = initialHeuristic;
		revAstarTestResult.m_algorithmInfo = "Rev-A*";
		TemplateAStar<state, action, environment> revAstar;
		state start = original;
		timer.StartTimer();
		bool solved = revAstar.GetPathTime(&env, goal, start, astarPath, SECONDS_LIMIT);
		timer.EndTimer();
		if (solved)
		{
			revAstarTestResult.m_solutionCost = env.GetPathLength(astarPath);
			revAstarTestResult.m_statesExpanded = revAstar.GetNodesExpanded();
			revAstarTestResult.m_neccessaryStatesExpanded = revAstar.GetNecessaryExpansions();
			revAstarTestResult.m_maxStatesInMemory = revAstar.getMemoryStatesUse();
		}
		revAstarTestResult.m_timeElapsed = timer.GetElapsedTime();
		testResults.push_back(revAstarTestResult);
	}
	// MM
	if (MMRun)
	{
		TestResult mmTestResult = TestResult(testInfo);
		mmTestResult.m_initialHeuristic = initialHeuristic;
		mmTestResult.m_algorithmInfo = "MM";
		MM<state, action, environment> mm;
		vector<state> mmPath;
		state start = original;
		timer.StartTimer();
		bool solved = mm.GetPath(&env, start, goal, &env, &env, mmPath, SECONDS_LIMIT);
		timer.EndTimer();
		if (solved)
		{
			mmTestResult.m_solutionCost = env.GetPathLength(mmPath);
			mmTestResult.m_statesExpanded = mm.GetNodesExpanded();
			mmTestResult.m_neccessaryStatesExpanded = mm.GetNecessaryExpansions();
			mmTestResult.m_maxStatesInMemory = mm.getMemoryStatesUse();
		}
		mmTestResult.m_timeElapsed = timer.GetElapsedTime();
		testResults.push_back(mmTestResult);
		if (printAbstractAlgos)
		{
			TestResult immTestResult = TestResult(testInfo);
			immTestResult.m_initialHeuristic = mmTestResult.m_initialHeuristic;
			immTestResult.m_algorithmInfo = "I-MM";
			immTestResult.m_solutionCost = mmTestResult.m_solutionCost;
			immTestResult.m_statesExpanded = mm.getIMMExpansions();
			immTestResult.m_timeElapsed = mmTestResult.m_timeElapsed;
			testResults.push_back(immTestResult);
		}
	}
	return testResults;
}

template <class state, class action, class environment, class hasher>
vector<TestResult> GenericTester<state, action, environment, hasher>::testLMA(state original, state goal, environment env, TestInfo testInfo)
{
	double initialHeuristic = env.HCost(original, goal);
	Timer timer;
	vector<state> idaPath;
	vector<TestResult> testResults;
	// IDA*
	if (IDAstarRun)
	{
		TestResult idaTestResult = TestResult(testInfo);
		idaTestResult.m_initialHeuristic = initialHeuristic;
		idaTestResult.m_algorithmInfo = "IDA*";
		IDAStar<state, action, false> idastar;
		state start = original;
		timer.StartTimer();
		bool solved = idastar.GetPath(&env, start, goal, idaPath, SECONDS_LIMIT);
		timer.EndTimer();
		if (solved)
		{
			idaTestResult.m_solutionCost = env.GetPathLength(idaPath);
			idaTestResult.m_statesExpanded = idastar.GetNodesExpanded();
			idaTestResult.m_neccessaryStatesExpanded = idastar.GetNecessaryExpansions();
		}
		idaTestResult.m_timeElapsed = timer.GetElapsedTime();
		testResults.push_back(idaTestResult);
		if (printAbstractAlgos)
		{
			TestResult dAstarTestResult = TestResult(testInfo);
			dAstarTestResult.m_initialHeuristic = idaTestResult.m_initialHeuristic;
			dAstarTestResult.m_algorithmInfo = "D-A*";
			dAstarTestResult.m_solutionCost = idaTestResult.m_solutionCost;
			dAstarTestResult.m_statesExpanded = idastar.getDAstarExpansions();
			dAstarTestResult.m_timeElapsed = idaTestResult.m_timeElapsed;
			testResults.push_back(dAstarTestResult);
		}
	}

	//IDBiHS
	if (IDBiHSRun)
	{
		TestResult testResult = TestResult(testInfo);
		testResult.m_initialHeuristic = initialHeuristic;
		testResult.m_algorithmInfo = "IDBiHS";
		IDBiHS<environment, state, action, false> idbihs(&env, F2Fheuristics, isConsistent, isUpdateByWorkload);
		state start = original;
		state midState;
		timer.StartTimer();
		bool solved = idbihs.GetMidState(start, goal, midState, SECONDS_LIMIT);
		timer.EndTimer();
		if (solved)
		{
			testResult.m_solutionCost = idbihs.getPathLength();
			testResult.m_statesExpanded = idbihs.GetNodesExpanded();
			testResult.m_neccessaryStatesExpanded = idbihs.GetNecessaryExpansions();
		}
		testResult.m_timeElapsed = timer.GetElapsedTime();
		testResults.push_back(testResult);

		if (printAbstractAlgos)
		{
			TestResult dMmtestResult = TestResult(testInfo);
			dMmtestResult.m_initialHeuristic = initialHeuristic;
			dMmtestResult.m_algorithmInfo = "D-MM";
			dMmtestResult.m_solutionCost = testResult.m_solutionCost;
			dMmtestResult.m_statesExpanded = idbihs.getDMMExpansions();
			dMmtestResult.m_timeElapsed = testResult.m_timeElapsed;
			testResults.push_back(dMmtestResult);
		}
	}
	return testResults;
}

template <class state, class action, class environment, class hasher>
vector<TestResult> GenericTester<state, action, environment, hasher>::testFMA(state original, state goal, environment env, TestInfo testInfo, unsigned long statesQuantityBound, vector<double> quantityPercentages)
{
	double initialHeuristic = env.HCost(original, goal);
	Timer timer;

	vector<state> astarPath;
	vector<state> idaPath;
	vector<TestResult> testResults;

	//FullBFBDS
	if (BFBDSRUN)
	{
		vector<state> fullBfbdsPath;
		state midState;
		bool solved;
		unsigned long nodesExpanded;
		for (double quantityPercentage : quantityPercentages)
		{
			TestResult bfbdsTestResult = TestResult(testInfo);
			bfbdsTestResult.m_initialHeuristic = initialHeuristic;
			bfbdsTestResult.m_algorithmInfo = "BFBDS-" + to_string(quantityPercentage) + "%";
			unsigned long statesQuantityBoundforBFBDS = std::max(statesQuantityBound * quantityPercentage, MINIMUM_STATES_QUANTITY_BOUND);
			state start = original;
			timer.StartTimer();
			BFBDS<state, action, environment, MbbdsBloomFilter<state, hasher>, false> bfbds(&env, statesQuantityBoundforBFBDS, isUpdateByWorkload, isConsistent, revAlgo, F2Fheuristics);
			solved = bfbds.GetMidState(start, goal, midState, fullBfbdsPath, SECONDS_LIMIT, threePhase);
			timer.EndTimer();

			bfbdsTestResult.m_maxStatesInMemory = statesQuantityBoundforBFBDS;
			if (solved)
			{
				bfbdsTestResult.m_solutionCost = bfbds.getPathLength();
				bfbdsTestResult.m_statesExpanded = bfbds.getNodesExpanded();
				bfbdsTestResult.m_neccessaryStatesExpanded = bfbds.getNecessaryExpansions();
				bfbdsTestResult.m_iterations = bfbds.getIterationsNum();
			}
			else
			{
				break;
			}
			bfbdsTestResult.m_timeElapsed = timer.GetElapsedTime();
			testResults.push_back(bfbdsTestResult);
		}
	}
	if (ASTARpIDBiHS)
	{
		for (double quantityPercentage : quantityPercentages)
		{
			TestResult testResult = TestResult(testInfo);
			testResult.m_initialHeuristic = initialHeuristic;
			testResult.m_algorithmInfo = "A*+IDBiHS-" + to_string(quantityPercentage) + "%";
			unsigned long statesQuantityBoundforASPIDBiHS = std::max(statesQuantityBound * quantityPercentage, MINIMUM_STATES_QUANTITY_BOUND);

			IDBiHS<environment, state, action, false> idbihs(&env, F2Fheuristics, isConsistent, isUpdateByWorkload);
			state midState;
			timer.StartTimer();
			state start = original;
			bool solved = idbihs.Astar_plus_IDBiHS(start, goal, midState, statesQuantityBoundforASPIDBiHS, SECONDS_LIMIT - timer.GetElapsedTime(), detectDuplicate);
			timer.EndTimer();
			testResult.m_maxStatesInMemory = statesQuantityBoundforASPIDBiHS;
			if (solved)
			{
				testResult.m_solutionCost = idbihs.getPathLength();
				testResult.m_statesExpanded = idbihs.GetNodesExpanded();
			}
			else
			{
				break;
			}
			testResult.m_timeElapsed = timer.GetElapsedTime();
			testResults.push_back(testResult);
		}
	}
	if (IDTHSpTrans)
	{
		{
			for (double quantityPercentage : quantityPercentages)
			{
				TestResult testResult = TestResult(testInfo);
				testResult.m_initialHeuristic = initialHeuristic;
				testResult.m_algorithmInfo = "IDTHSpTrans-" + to_string(quantityPercentage) + "%";
				unsigned long statesQuantityBoundforASPIDBiHS = std::max(statesQuantityBound * quantityPercentage, MINIMUM_STATES_QUANTITY_BOUND);
				state start = original;
				timer.StartTimer();
				unsigned long nodesExpanded = 0;
				unsigned long necessaryNodesExpanded = 0;
				IDTHSwTrans<state, action, false> idbihs(F2Fheuristics, isConsistent, isUpdateByWorkload, 1, true);
				state midState;
				bool solved = idbihs.GetPath(&env, start, goal, SECONDS_LIMIT, statesQuantityBoundforASPIDBiHS);
				nodesExpanded += idbihs.GetNodesExpanded();
				necessaryNodesExpanded += idbihs.GetNecessaryExpansions();
				timer.EndTimer();
				testResult.m_maxStatesInMemory = statesQuantityBoundforASPIDBiHS;
				if (solved)
				{
					testResult.m_solutionCost = idbihs.getPathLength();
					testResult.m_statesExpanded = nodesExpanded;
					testResult.m_neccessaryStatesExpanded = necessaryNodesExpanded;
				}
				else
				{
					break;
				}
				testResult.m_timeElapsed = timer.GetElapsedTime();
				testResults.push_back(testResult);
			}
		}
		{
			for (double quantityPercentage : quantityPercentages)
			{
				TestResult testResult = TestResult(testInfo);
				testResult.m_initialHeuristic = initialHeuristic;
				testResult.m_algorithmInfo = "IDTHSpTrans_NDD-" + to_string(quantityPercentage) + "%";
				unsigned long statesQuantityBoundforASPIDBiHS = std::max(statesQuantityBound * quantityPercentage, MINIMUM_STATES_QUANTITY_BOUND);
				state start = original;
				timer.StartTimer();
				unsigned long nodesExpanded = 0;
				unsigned long necessaryNodesExpanded = 0;
				IDTHSwTrans<state, action, false> idbihs(F2Fheuristics, isConsistent, isUpdateByWorkload, 1, false);
				state midState;
				bool solved = idbihs.GetPath(&env, start, goal, SECONDS_LIMIT, statesQuantityBoundforASPIDBiHS);
				nodesExpanded += idbihs.GetNodesExpanded();
				necessaryNodesExpanded += idbihs.GetNecessaryExpansions();
				timer.EndTimer();
				testResult.m_maxStatesInMemory = statesQuantityBoundforASPIDBiHS;
				if (solved)
				{
					testResult.m_solutionCost = idbihs.getPathLength();
					testResult.m_statesExpanded = nodesExpanded;
					testResult.m_neccessaryStatesExpanded = necessaryNodesExpanded;
				}
				else
				{
					break;
				}
				testResult.m_timeElapsed = timer.GetElapsedTime();
				testResults.push_back(testResult);
			}
		}
	}
	if (AstarPIDAstarRun)
	{
		for (double quantityPercentage : quantityPercentages)
		{
			TestResult testResult = TestResult(testInfo);
			testResult.m_initialHeuristic = initialHeuristic;
			testResult.m_algorithmInfo = "Astar+IDAstar-" + to_string(quantityPercentage) + "%";
			unsigned long statesQuantityBoundforASPIDAS = std::max(statesQuantityBound * quantityPercentage, MINIMUM_STATES_QUANTITY_BOUND);
			TemplateAStar<state, action, environment> astar;
			state start = original;
			timer.StartTimer();
			bool solved = astar.GetPathTime(&env, start, goal, astarPath, SECONDS_LIMIT, true, statesQuantityBoundforASPIDAS);
			unsigned long nodesExpanded = astar.GetNodesExpanded();
			unsigned long necessaryNodesExpanded = 0;
			testResult.m_maxStatesInMemory = statesQuantityBoundforASPIDAS;
			if (solved)
			{
				timer.EndTimer();
				necessaryNodesExpanded = astar.GetNecessaryExpansions();
				testResult.m_solutionCost = env.GetPathLength(astarPath);
				testResult.m_statesExpanded = nodesExpanded;
				testResult.m_neccessaryStatesExpanded = necessaryNodesExpanded;
			}
			else
			{
				IDAStar<state, action, false> idastar;
				solved = idastar.ASpIDA(&env, start, goal, idaPath, astar.getStatesList(), SECONDS_LIMIT - timer.GetElapsedTime(), detectDuplicate);
				nodesExpanded += idastar.GetNodesExpanded();
				necessaryNodesExpanded += idastar.GetNecessaryExpansions();
				timer.EndTimer();
				if (solved)
				{
					necessaryNodesExpanded = astar.GetNecessaryExpansions();
					testResult.m_solutionCost = idastar.getSolLength();
					testResult.m_statesExpanded = nodesExpanded;
					testResult.m_neccessaryStatesExpanded = necessaryNodesExpanded;
				}
				else
				{
					break;
				}
			}
			testResult.m_timeElapsed = timer.GetElapsedTime();
			testResults.push_back(testResult);
		}
	}
	if (AstarPIDAstarReverseRun)
	{
		for (double quantityPercentage : quantityPercentages)
		{
			TestResult testResult = TestResult(testInfo);
			testResult.m_initialHeuristic = initialHeuristic;
			testResult.m_algorithmInfo = "Astar+IDAstar+Reverse-" + to_string(quantityPercentage) + "%";
			unsigned long statesQuantityBoundforASPIDARS = std::max(statesQuantityBound * quantityPercentage, MINIMUM_STATES_QUANTITY_BOUND);
			TemplateAStar<state, action, environment> astar;
			state start = original;
			timer.StartTimer();
			bool solved = astar.GetPathTime(&env, goal, start, astarPath, SECONDS_LIMIT, true, statesQuantityBoundforASPIDARS, false);
			unsigned long nodesExpanded = astar.GetNodesExpanded();
			unsigned long necessaryNodesExpanded = 0;
			testResult.m_maxStatesInMemory = statesQuantityBoundforASPIDARS;
			if (solved)
			{
				timer.EndTimer();
				necessaryNodesExpanded = astar.GetNecessaryExpansions();
				testResult.m_solutionCost = env.GetPathLength(astarPath);
				testResult.m_statesExpanded = nodesExpanded;
				testResult.m_neccessaryStatesExpanded = necessaryNodesExpanded;
			}
			else
			{
				IDAStar<state, action, false> idastar;
				solved = idastar.ASpIDArev(&env, start, goal, idaPath, astar.getStatesList(), astar.getPrevF(), SECONDS_LIMIT - timer.GetElapsedTime());
				nodesExpanded += idastar.GetNodesExpanded();
				necessaryNodesExpanded += idastar.GetNecessaryExpansions();
				timer.EndTimer();
				if (solved)
				{
					testResult.m_solutionCost = idastar.getSolLength();
					testResult.m_statesExpanded = nodesExpanded;
					testResult.m_neccessaryStatesExpanded = necessaryNodesExpanded;
				}
				else
				{
					break;
				}
			}
			testResult.m_timeElapsed = timer.GetElapsedTime();
			testResults.push_back(testResult);
		}
	}
	if (AstarPIDAstarReverseMinHRun)
	{
		for (double quantityPercentage : quantityPercentages)
		{
			TestResult testResult = TestResult(testInfo);
			testResult.m_initialHeuristic = initialHeuristic;
			testResult.m_algorithmInfo = "Astar+IDAstar+Reverse-" + to_string(quantityPercentage) + "%";
			unsigned long statesQuantityBoundforASPIDARS = std::max(statesQuantityBound * quantityPercentage, MINIMUM_STATES_QUANTITY_BOUND);
			TemplateAStar<state, action, environment> astar;
			state start = original;
			timer.StartTimer();
			bool solved = astar.GetPathTime(&env, goal, start, astarPath, SECONDS_LIMIT, true, statesQuantityBoundforASPIDARS, false);
			unsigned long nodesExpanded = astar.GetNodesExpanded();
			unsigned long necessaryNodesExpanded = 0;
			testResult.m_maxStatesInMemory = statesQuantityBoundforASPIDARS;
			if (solved)
			{
				timer.EndTimer();
				necessaryNodesExpanded = astar.GetNecessaryExpansions();
				testResult.m_solutionCost = env.GetPathLength(astarPath);
				testResult.m_statesExpanded = nodesExpanded;
				testResult.m_neccessaryStatesExpanded = necessaryNodesExpanded;
			}
			else
			{
				IDAStar<state, action, false> idastar;
				solved = idastar.ASpIDArev(&env, start, goal, idaPath, astar.getStatesList(), astar.getPrevF(), SECONDS_LIMIT - timer.GetElapsedTime(), isConsistent, true);
				nodesExpanded += idastar.GetNodesExpanded();
				necessaryNodesExpanded += idastar.GetNecessaryExpansions();
				timer.EndTimer();
				if (solved)
				{
					testResult.m_solutionCost = idastar.getSolLength();
					testResult.m_statesExpanded = nodesExpanded;
					testResult.m_neccessaryStatesExpanded = necessaryNodesExpanded;
				}
				else
				{
					break;
				}
			}
			testResult.m_timeElapsed = timer.GetElapsedTime();
			testResults.push_back(testResult);
		}
	}
	if (BAI)
	{
		for (double quantityPercentage : quantityPercentages)
		{
			TestResult testResult = TestResult(testInfo);
			testResult.m_initialHeuristic = initialHeuristic;
			testResult.m_algorithmInfo = "BAI-" + to_string(quantityPercentage) + "%";
			unsigned long statesQuantityBoundforASPIDARS = std::max(statesQuantityBound * quantityPercentage, MINIMUM_STATES_QUANTITY_BOUND);
			TemplateAStar<state, action, environment> astar;
			state start = original;
			timer.StartTimer();
			bool solved = astar.GetPathTime(&env, goal, start, astarPath, SECONDS_LIMIT, true, statesQuantityBoundforASPIDARS);
			unsigned long nodesExpanded = astar.GetNodesExpanded();
			unsigned long necessaryNodesExpanded = 0;
			testResult.m_maxStatesInMemory = statesQuantityBoundforASPIDARS;
			if (solved)
			{
				timer.EndTimer();
				necessaryNodesExpanded = astar.GetNecessaryExpansions();
				testResult.m_solutionCost = env.GetPathLength(astarPath);
				testResult.m_statesExpanded = nodesExpanded;
				testResult.m_neccessaryStatesExpanded = necessaryNodesExpanded;
			}
			else
			{
				IDAStar<state, action, false> idastar;
				solved = idastar.BAI(&env, start, goal, idaPath, astar.getStatesList(), SECONDS_LIMIT - timer.GetElapsedTime(), false);
				nodesExpanded += idastar.GetNodesExpanded();
				necessaryNodesExpanded += idastar.GetNecessaryExpansions();
				timer.EndTimer();
				if (solved)
				{
					testResult.m_solutionCost = idastar.getSolLength();
					testResult.m_statesExpanded = nodesExpanded;
					testResult.m_neccessaryStatesExpanded = necessaryNodesExpanded;
				}
				else
				{
					break;
				}
			}
			testResult.m_timeElapsed = timer.GetElapsedTime();
			testResults.push_back(testResult);
		}
	}
	if (Max_BAI)
	{
		for (double quantityPercentage : quantityPercentages)
		{
			TestResult testResult = TestResult(testInfo);
			testResult.m_initialHeuristic = initialHeuristic;
			testResult.m_algorithmInfo = "MaxBAI-" + to_string(quantityPercentage) + "%";
			unsigned long statesQuantityBoundforASPIDARS = std::max(statesQuantityBound * quantityPercentage, MINIMUM_STATES_QUANTITY_BOUND);
			TemplateAStar<state, action, environment> astar;
			state start = original;
			timer.StartTimer();
			bool solved = astar.GetPathTime(&env, goal, start, astarPath, SECONDS_LIMIT, true, statesQuantityBoundforASPIDARS);
			unsigned long nodesExpanded = astar.GetNodesExpanded();
			unsigned long necessaryNodesExpanded = 0;
			testResult.m_maxStatesInMemory = statesQuantityBoundforASPIDARS;
			if (solved)
			{
				timer.EndTimer();
				necessaryNodesExpanded = astar.GetNecessaryExpansions();
				testResult.m_solutionCost = env.GetPathLength(astarPath);
				testResult.m_statesExpanded = nodesExpanded;
				testResult.m_neccessaryStatesExpanded = necessaryNodesExpanded;
			}
			else
			{
				IDAStar<state, action, false> idastar;
				solved = idastar.BAI(&env, start, goal, idaPath, astar.getStatesList(), SECONDS_LIMIT - timer.GetElapsedTime(), true);
				nodesExpanded += idastar.GetNodesExpanded();
				necessaryNodesExpanded += idastar.GetNecessaryExpansions();
				timer.EndTimer();
				if (solved)
				{
					testResult.m_solutionCost = idastar.getSolLength();
					testResult.m_statesExpanded = nodesExpanded;
					testResult.m_neccessaryStatesExpanded = necessaryNodesExpanded;
				}
				else
				{
					break;
				}
			}
			testResult.m_timeElapsed = timer.GetElapsedTime();
			testResults.push_back(testResult);
		}
		return testResults;
	}
}

#endif
