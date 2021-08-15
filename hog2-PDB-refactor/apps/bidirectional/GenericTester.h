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
#include "GenericBloomFilter.h"
#include "TestInfo.h"
#include "TestResult.h"
#include "json.h"
#include <fstream>
#include <iostream>
#include <boost/format.hpp>
#include <math.h>
#include <ctime>

using namespace std;
using json = nlohmann::json;

template <class state, class action, class environment>
class GenericTester
{
public:
	GenericTester() {}
	virtual ~GenericTester() {}
	void genericTest(state original, state goal, environment env, ofstream &myfile, TestInfo testInfo, double smallestEdge, json configurations);
	vector<TestResult> testUMA(state original, state goal, environment env, TestInfo testInfo, json configurations);
	vector<TestResult> testLMA(state original, state goal, environment env, TestInfo testInfo, json configurations);
	vector<TestResult> testFMA(state original, state goal, environment env, TestInfo testInfo, unsigned long statesQuantityBound, double smallestEdge, json configurations);
};

template <class state, class action, class environment>
void GenericTester<state, action, environment>::genericTest(state original, state goal, environment env, ofstream &myfile, TestInfo testInfo, double smallestEdge, json configurations)
{
	vector<TestResult> lmaResults = testLMA(original, goal, env, testInfo, configurations);

	vector<TestResult> umaResults = testUMA(original, goal, env, testInfo, configurations);
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
		statesQuantityBound = configurations["properties"]["defaultStatesQuantityBound"].get<unsigned long>();
	}
	vector<TestResult> fmaResults = testFMA(original, goal, env, testInfo, statesQuantityBound, smallestEdge, configurations);

	vector<vector<TestResult>> testResultsVectorts{lmaResults, umaResults, fmaResults};
	for (vector<TestResult> testResultsVector : testResultsVectorts)
	{
		for (TestResult testResult : testResultsVector)
		{
			myfile << testResult.csvSerialize() << std::endl;
		}
	}
}

template <class state, class action, class environment>
vector<TestResult> GenericTester<state, action, environment>::testUMA(state original, state goal, environment env, TestInfo testInfo, json configurations)
{
	double initialHeuristic = env.HCost(original, goal);
	Timer timer;
	vector<state> astarPath;

	vector<TestResult> testResults;
	// A*
	vector<AStarOpenClosedDataWithF<state>> astarOpenList;
	if (configurations["algorithms"]["astar"].get<bool>())
	{
		TestResult astarTestResult = TestResult(testInfo);
		astarTestResult.m_initialHeuristic = initialHeuristic;
		astarTestResult.m_algorithmInfo = "A*";

		TemplateAStar<state, action, environment> astar;
		state start = original;
		timer.StartTimer();
		bool solved = astar.GetPathTime(&env, start, goal, astarPath, configurations["properties"]["secondsLimit"].get<int>());
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
		if (configurations["properties"]["printAbstractAlgos"].get<bool>())
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
	if (configurations["algorithms"]["reversedAstar"].get<bool>())
	{
		TestResult revAstarTestResult = TestResult(testInfo);
		revAstarTestResult.m_initialHeuristic = initialHeuristic;
		revAstarTestResult.m_algorithmInfo = "Rev-A*";
		TemplateAStar<state, action, environment> revAstar;
		state start = original;
		timer.StartTimer();
		bool solved = revAstar.GetPathTime(&env, goal, start, astarPath, configurations["properties"]["secondsLimit"].get<int>());
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
	if (configurations["algorithms"]["mm"].get<bool>())
	{
		TestResult mmTestResult = TestResult(testInfo);
		mmTestResult.m_initialHeuristic = initialHeuristic;
		mmTestResult.m_algorithmInfo = "MM";
		MM<state, action, environment> mm;
		vector<state> mmPath;
		state start = original;
		timer.StartTimer();
		bool solved = mm.GetPath(&env, start, goal, &env, &env, mmPath, configurations["properties"]["secondsLimit"].get<int>());
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
		if (configurations["properties"]["printAbstractAlgos"].get<bool>())
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

template <class state, class action, class environment>
vector<TestResult> GenericTester<state, action, environment>::testLMA(state original, state goal, environment env, TestInfo testInfo, json configurations)
{
	double initialHeuristic = env.HCost(original, goal);
	Timer timer;
	vector<state> idaPath;
	vector<TestResult> testResults;
	// IDA*
	if (configurations["algorithms"]["idastar"].get<bool>())
	{
		TestResult idaTestResult = TestResult(testInfo);
		idaTestResult.m_initialHeuristic = initialHeuristic;
		idaTestResult.m_algorithmInfo = "IDA*";
		IDAStar<state, action, false> idastar;
		state start = original;
		timer.StartTimer();
		bool solved = idastar.GetPath(&env, start, goal, idaPath, configurations["properties"]["secondsLimit"].get<int>());
		timer.EndTimer();
		if (solved)
		{
			idaTestResult.m_solutionCost = env.GetPathLength(idaPath);
			idaTestResult.m_statesExpanded = idastar.GetNodesExpanded();
			idaTestResult.m_neccessaryStatesExpanded = idastar.GetNecessaryExpansions();
		}
		idaTestResult.m_timeElapsed = timer.GetElapsedTime();
		testResults.push_back(idaTestResult);
		if (configurations["properties"]["printAbstractAlgos"].get<bool>())
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
	if (configurations["algorithms"]["idbihs"].get<bool>())
	{
		TestResult testResult = TestResult(testInfo);
		testResult.m_initialHeuristic = initialHeuristic;
		testResult.m_algorithmInfo = "IDBiHS";
		IDBiHS<environment, state, action, false> idbihs(&env, configurations["properties"]["f2fHeuristics"].get<bool>(), configurations["properties"]["consistent"].get<bool>(), configurations["properties"]["UpdateByWorkload"].get<bool>());
		state start = original;
		state midState;
		timer.StartTimer();
		bool solved = idbihs.GetMidState(start, goal, midState, configurations["properties"]["secondsLimit"].get<int>());
		timer.EndTimer();
		if (solved)
		{
			testResult.m_solutionCost = idbihs.getPathLength();
			testResult.m_statesExpanded = idbihs.GetNodesExpanded();
			testResult.m_neccessaryStatesExpanded = idbihs.GetNecessaryExpansions();
		}
		testResult.m_timeElapsed = timer.GetElapsedTime();
		testResults.push_back(testResult);

		if (configurations["properties"]["printAbstractAlgos"].get<bool>())
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

template <class state, class action, class environment>
vector<TestResult> GenericTester<state, action, environment>::testFMA(state original, state goal, environment env, TestInfo testInfo, unsigned long statesQuantityBound, double smallestEdge, json configurations)
{
	vector<double> quantityPercentages = configurations["properties"]["rmaStatesPercentages"].get<vector<double>>();
	double initialHeuristic = env.HCost(original, goal);
	Timer timer;

	vector<state> astarPath;
	vector<state> idaPath;
	vector<TestResult> testResults;

	//FullBFBDS
	if (configurations["algorithms"]["bfbds"].get<bool>())
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
			unsigned long statesQuantityBoundforBFBDS = std::max(statesQuantityBound * quantityPercentage, configurations["properties"]["minStatesQuantityBound"].get<double>());
			state start = original;
			timer.StartTimer();
			BFBDS<state, action, environment> bfbds(&env, statesQuantityBoundforBFBDS, configurations["properties"]["UpdateByWorkload"].get<bool>(), configurations["properties"]["consistent"].get<bool>(), configurations["properties"]["bfbdsReverse"].get<bool>(), configurations["properties"]["f2fHeuristics"].get<bool>(), configurations["properties"]["bfbdsVerbose"].get<bool>(), smallestEdge);
			solved = bfbds.GetMidState(start, goal, midState, fullBfbdsPath, configurations["properties"]["secondsLimit"].get<int>(), configurations["properties"]["bfbfsThreePhase"].get<bool>());
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
	if (configurations["algorithms"]["astarNidbihs"].get<bool>())
	{
		for (double quantityPercentage : quantityPercentages)
		{
			TestResult testResult = TestResult(testInfo);
			testResult.m_initialHeuristic = initialHeuristic;
			testResult.m_algorithmInfo = "A*+IDBiHS-" + to_string(quantityPercentage) + "%";
			unsigned long statesQuantityBoundforASPIDBiHS = std::max(statesQuantityBound * quantityPercentage, configurations["properties"]["minStatesQuantityBound"].get<double>());

			IDBiHS<environment, state, action, false> idbihs(&env, configurations["properties"]["f2fHeuristics"].get<bool>(), configurations["properties"]["consistent"].get<bool>(), configurations["properties"]["UpdateByWorkload"].get<bool>());
			state midState;
			timer.StartTimer();
			state start = original;
			bool solved = idbihs.Astar_plus_IDBiHS(start, goal, midState, statesQuantityBoundforASPIDBiHS, configurations["properties"]["secondsLimit"].get<int>() - timer.GetElapsedTime(), configurations["properties"]["detectDuplicates"].get<bool>());
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
	if (configurations["algorithms"]["idthsTrans"].get<bool>())
	{
		{
			for (double quantityPercentage : quantityPercentages)
			{
				TestResult testResult = TestResult(testInfo);
				testResult.m_initialHeuristic = initialHeuristic;
				testResult.m_algorithmInfo = "idthsTrans-" + to_string(quantityPercentage) + "%";
				unsigned long statesQuantityBoundforASPIDBiHS = std::max(statesQuantityBound * quantityPercentage, configurations["properties"]["minStatesQuantityBound"].get<double>());
				state start = original;
				timer.StartTimer();
				unsigned long nodesExpanded = 0;
				unsigned long necessaryNodesExpanded = 0;
				IDTHSwTrans<state, action, false> idbihs(configurations["properties"]["f2fHeuristics"].get<bool>(), configurations["properties"]["consistent"].get<bool>(), configurations["properties"]["UpdateByWorkload"].get<bool>(), 1, false);
				state midState;
				bool solved = idbihs.GetPath(&env, start, goal, configurations["properties"]["secondsLimit"].get<int>(), statesQuantityBoundforASPIDBiHS);
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
				unsigned long statesQuantityBoundforASPIDBiHS = std::max(statesQuantityBound * quantityPercentage, configurations["properties"]["minStatesQuantityBound"].get<double>());
				state start = original;
				timer.StartTimer();
				unsigned long nodesExpanded = 0;
				unsigned long necessaryNodesExpanded = 0;
				IDTHSwTrans<state, action, false> idbihs(configurations["properties"]["f2fHeuristics"].get<bool>(), configurations["properties"]["consistent"].get<bool>(), configurations["properties"]["UpdateByWorkload"].get<bool>(), 1, false);
				state midState;
				bool solved = idbihs.GetPath(&env, start, goal, configurations["properties"]["secondsLimit"].get<int>(), statesQuantityBoundforASPIDBiHS);
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
	if (configurations["algorithms"]["astarNidastar"].get<bool>())
	{
		for (double quantityPercentage : quantityPercentages)
		{
			TestResult testResult = TestResult(testInfo);
			testResult.m_initialHeuristic = initialHeuristic;
			testResult.m_algorithmInfo = "Astar+IDAstar-" + to_string(quantityPercentage) + "%";
			unsigned long statesQuantityBoundforASPIDAS = std::max(statesQuantityBound * quantityPercentage, configurations["properties"]["minStatesQuantityBound"].get<double>());
			TemplateAStar<state, action, environment> astar;
			state start = original;
			timer.StartTimer();
			bool solved = astar.GetPathTime(&env, start, goal, astarPath, configurations["properties"]["secondsLimit"].get<int>(), true, statesQuantityBoundforASPIDAS);
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
				solved = idastar.ASpIDA(&env, start, goal, idaPath, astar.getStatesList(), configurations["properties"]["secondsLimit"].get<int>() - timer.GetElapsedTime(), configurations["properties"]["detectDuplicates"].get<bool>());
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
	if (configurations["algorithms"]["reversedAstarNidastar"].get<bool>())
	{
		for (double quantityPercentage : quantityPercentages)
		{
			TestResult testResult = TestResult(testInfo);
			testResult.m_initialHeuristic = initialHeuristic;
			testResult.m_algorithmInfo = "Astar+IDAstar+Reverse-" + to_string(quantityPercentage) + "%";
			unsigned long statesQuantityBoundforASPIDARS = std::max(statesQuantityBound * quantityPercentage, configurations["properties"]["minStatesQuantityBound"].get<double>());
			TemplateAStar<state, action, environment> astar;
			state start = original;
			timer.StartTimer();
			bool solved = astar.GetPathTime(&env, goal, start, astarPath, configurations["properties"]["secondsLimit"].get<int>(), true, statesQuantityBoundforASPIDARS, false);
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
				solved = idastar.ASpIDArev(&env, start, goal, idaPath, astar.getStatesList(), astar.getPrevF(), configurations["properties"]["secondsLimit"].get<int>() - timer.GetElapsedTime());
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
	if (configurations["algorithms"]["reversedAstarNidastarMinHRun"].get<bool>())
	{
		for (double quantityPercentage : quantityPercentages)
		{
			TestResult testResult = TestResult(testInfo);
			testResult.m_initialHeuristic = initialHeuristic;
			testResult.m_algorithmInfo = "Astar+IDAstar+ReverseMinH-" + to_string(quantityPercentage) + "%";
			unsigned long statesQuantityBoundforASPIDARS = std::max(statesQuantityBound * quantityPercentage, configurations["properties"]["minStatesQuantityBound"].get<double>());
			TemplateAStar<state, action, environment> astar;
			state start = original;
			timer.StartTimer();
			bool solved = astar.GetPathTime(&env, goal, start, astarPath, configurations["properties"]["secondsLimit"].get<int>(), true, statesQuantityBoundforASPIDARS, false);
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
				solved = idastar.ASpIDArev(&env, start, goal, idaPath, astar.getStatesList(), astar.getPrevF(), configurations["properties"]["secondsLimit"].get<int>() - timer.GetElapsedTime(), configurations["properties"]["consistent"].get<bool>(), true);
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
	if (configurations["algorithms"]["bai"].get<bool>())
	{
		for (double quantityPercentage : quantityPercentages)
		{
			TestResult testResult = TestResult(testInfo);
			testResult.m_initialHeuristic = initialHeuristic;
			testResult.m_algorithmInfo = "BAI-" + to_string(quantityPercentage) + "%";
			unsigned long statesQuantityBoundforASPIDARS = std::max(statesQuantityBound * quantityPercentage, configurations["properties"]["minStatesQuantityBound"].get<double>());
			TemplateAStar<state, action, environment> astar;
			state start = original;
			timer.StartTimer();
			bool solved = astar.GetPathTime(&env, goal, start, astarPath, configurations["properties"]["secondsLimit"].get<int>(), true, statesQuantityBoundforASPIDARS);
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
				solved = idastar.BAI(&env, start, goal, idaPath, astar.getStatesList(), configurations["properties"]["secondsLimit"].get<int>() - timer.GetElapsedTime(), false);
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
	if (configurations["algorithms"]["maxBai"].get<bool>())
	{
		for (double quantityPercentage : quantityPercentages)
		{
			TestResult testResult = TestResult(testInfo);
			testResult.m_initialHeuristic = initialHeuristic;
			testResult.m_algorithmInfo = "MaxBAI-" + to_string(quantityPercentage) + "%";
			unsigned long statesQuantityBoundforASPIDARS = std::max(statesQuantityBound * quantityPercentage, configurations["properties"]["minStatesQuantityBound"].get<double>());
			TemplateAStar<state, action, environment> astar;
			state start = original;
			timer.StartTimer();
			bool solved = astar.GetPathTime(&env, goal, start, astarPath, configurations["properties"]["secondsLimit"].get<int>(), true, statesQuantityBoundforASPIDARS);
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
				solved = idastar.BAI(&env, start, goal, idaPath, astar.getStatesList(), configurations["properties"]["secondsLimit"].get<int>() - timer.GetElapsedTime(), true);
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
