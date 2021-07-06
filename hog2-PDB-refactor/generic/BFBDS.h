/*
 *  BFBDS.h
 *
 *  Created by Steven Danishevski on DEC 19.
 *
 */

#ifndef BFBDS_H
#define BFBDS_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <math.h>
#include "MM.h"
#include "IDBiHS.h"
#include "Bloom.h"

const bool DEFAULT_UPDATE_BY_WORKLOAD = true;
const bool DEFAULT_IS_CONSISTENT = false;
const bool DEFAULT_REV_ALGO = false;
const bool DEFAULT_F2F_HEURISTICS = false;
const bool DEFAULT_F2F_VERBOSE = false;
const bool DEFAULT_SECONDS_LIMIT = 60 * 10;
const double DEFAULT_SMALLEST_EDGE = 1;
const int SATURATION_MAX_INC = 1;
const int DEFAULT_HASH_NUM = 5;
const int MIN_HASH_NUM = 1;
const int MAX_HASH_NUM = 10;
const double LN2 = 0.693;

template <class state, class action, class environment>
class BFBDS
{
	typedef std::function<bool(const state &, const state &)> StatesComparator;

public:
	BFBDS(environment *env, unsigned long statesQuantityBound, bool isUpdateByWorkload = DEFAULT_UPDATE_BY_WORKLOAD, bool isConsistent = DEFAULT_IS_CONSISTENT, bool revAlgo = DEFAULT_REV_ALGO, bool F2Fheuristics = DEFAULT_F2F_HEURISTICS, bool DEFAULT_F2F_VERBOSE = false, double smallestEdge = DEFAULT_SMALLEST_EDGE);
	virtual ~BFBDS() { deleteBloomFilters(); }
	bool GetMidState(state originState, state goalState, state &midState, std::vector<state> &thePath, int secondsLimit = 600, bool threePhase = true);
	double getPathLength() { return pathLength; }
	uint64_t getNodesExpanded() { return nodesExpanded; }
	uint64_t getNecessaryExpansions() { return necessaryExpansions; }
	int getIterationsNum() { return iteration_num; }

private:
	bool GetMidState(state originState, state goalState, state &midState, int secondsLimit = DEFAULT_SECONDS_LIMIT, double startingFBound = 0);
	bool DoIteration(state &originState, state &goalState, state &parent, state &currState, double currentDirectionBound, double g, state &midState);
	bool checkState(state &midState, double g);
	void updateBounds();
	double calculateNextForwardBound();
	double isNextForwardSearch();
	void UpdateNextBound(double potentialNextBound);
	double calculateBoundByWorkload(double newBound, double previousBound, double currentDirectionBound, uint64_t currentLoad, uint64_t currentOppositeLoad);
	int calculateOptimalHashNum(uint64_t bloomfilterSize, uint64_t expectedEntries);
	void resetBfSearchValues();
	void deleteBloomFilters();
	void printIterationInfo();
	BloomFilter *getNextBloomfilter();

	//Constructor variables
	environment *env;
	bool threePhase;
	bool isUpdateByWorkload;
	bool isConsistent;
	bool revAlgo;
	bool F2Fheuristics;
	uint64_t statesQuantityBound;
	double smallestEdge;

	//Results variables
	uint64_t nodesExpanded, nodesTouched, necessaryExpansions;
	int iteration_num;
	double pathLength;

	//Run variables
	bool verbose;
	unsigned long memoryBound;
	BloomFilter *previousBloomfilter = nullptr;
	BloomFilter *currentBloomfilter = nullptr;
	std::map<state, double, StatesComparator> middleStates;
	StatesComparator statesComparator;
	bool forwardSearch;
	bool listReady;
	bool outOfSpace;
	bool firstRun;
	double backwardBound, forwardBound, fBound;
	double previousMaxG;
	double currentMaxG;
	unsigned long lastIterBoundExpansions;
	unsigned long forwardExpandedInLastIter;
	unsigned long backwardExpandedInLastIter;
	double minCurrentError = std::numeric_limits<double>::max();
	double minPreviousError = 0;

	//Next Run Variables
	double nextBound;
};

template <class state, class action, class environment>
BFBDS<state, action, environment>::BFBDS(environment *env, unsigned long statesQuantityBound, bool isUpdateByWorkload, bool isConsistent, bool revAlgo, bool F2Fheuristics, bool verbose, double smallestEdge) : env(env), statesQuantityBound(statesQuantityBound), isUpdateByWorkload(isUpdateByWorkload), isConsistent(isConsistent), revAlgo(revAlgo), F2Fheuristics(F2Fheuristics), verbose(verbose), smallestEdge(smallestEdge)
{
	this->statesComparator = [&](const state &firstState, const state &secondState)
	{ return env->GetStateHash(firstState) < env->GetStateHash(secondState); };
	this->middleStates = std::map<state, double, StatesComparator>(this->statesComparator);
}

template <class state, class action, class environment>
bool BFBDS<state, action, environment>::GetMidState(state originState, state goalState, state &midState, std::vector<state> &thePath, int secondsLimit, bool threePhase)
{
	if (this->verbose)
	{
		printf("Starting to solve with BFBDS, memory limit= %1.1llu\n", this->memoryBound);
	}
	this->nodesExpanded = this->nodesTouched = this->necessaryExpansions = 0;
	double lastBound;
	auto startTime = std::chrono::steady_clock::now();
	double elapsed_seconds;
	if (threePhase)
	{
		if (this->verbose)
		{
			printf("Starting BFBDS: MM phase\n");
		}
		MM<state, action, environment> mm;
		elapsed_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();
		bool solved = mm.GetPath(this->env, originState, goalState, this->env, this->env, thePath, secondsLimit - elapsed_seconds, statesQuantityBound);
		this->nodesExpanded = mm.GetNodesExpanded();
		this->necessaryExpansions = mm.GetNecessaryExpansions();
		lastBound = mm.getLastBound();
		if (solved)
		{
			this->pathLength = this->env->GetPathLength(thePath);
			return true;
		}
	}
	if (this->verbose)
	{
		printf("Starting BFBDS: BF phase\n");
	}
	elapsed_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();
	bool solved = GetMidState(originState, goalState, midState, secondsLimit - elapsed_seconds, lastBound);
	deleteBloomFilters();
	lastBound = this->fBound;
	if (!solved)
	{
		if (this->verbose)
		{
			printf("Starting BFBDS: IDBiHS phase\n");
		}
		IDBiHS<environment, state, action, false> idbihs(this->env, this->F2Fheuristics, this->isConsistent, this->isUpdateByWorkload, this->smallestEdge);
		elapsed_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();
		solved = idbihs.GetMidState(originState, goalState, midState, secondsLimit - elapsed_seconds, lastBound);
		unsigned long test = this->nodesExpanded;
		this->nodesExpanded += idbihs.GetNodesExpanded();
		if (solved)
		{
			this->pathLength = idbihs.getPathLength();
			this->necessaryExpansions = idbihs.GetNecessaryExpansions();
		}
	}
	return solved;
}

template <class state, class action, class environment>
double BFBDS<state, action, environment>::calculateBoundByWorkload(double newBound, double previousBound, double currentDirectionBound, uint64_t currentLoad, uint64_t currentOppositeLoad)
{
	if (currentLoad <= currentOppositeLoad)
	{
		return currentDirectionBound + newBound - previousBound;
	}
	else
	{
		return currentDirectionBound;
	}
}

template <class state, class action, class environment>
bool BFBDS<state, action, environment>::GetMidState(state originState, state goalState, state &midState, int secondsLimit, double startingFBound)
{
	state currentOriginState, currentGoalState;
	this->memoryBound = sizeof(originState) * statesQuantityBound;
	double initialHeuristic = this->env->HCost(originState, goalState);
	this->nextBound = this->fBound = std::max(smallestEdge, std::max(initialHeuristic, startingFBound));
	this->forwardSearch = true;
	this->forwardBound = this->fBound / 2;
	this->backwardBound = this->fBound - this->forwardBound;
	this->iteration_num = 0;
	this->lastIterBoundExpansions = 0;
	auto startTime = std::chrono::steady_clock::now();
	while (true)
	{
		resetBfSearchValues();
		int saturationIncreasedCount = 0;
		double last_saturation;
		while (true)
		{
			this->previousBloomfilter = this->currentBloomfilter;
			this->currentBloomfilter = nullptr;
			this->previousMaxG = this->currentMaxG;
			this->currentMaxG = 0;

			auto currentTime = std::chrono::steady_clock::now();
			std::chrono::duration<double> elapsed_seconds = currentTime - startTime;
			if (elapsed_seconds.count() >= secondsLimit)
			{
				return false;
			}

			if (this->verbose)
			{
				printIterationInfo();
			}

			//Initializing iteration variables
			this->outOfSpace = false;
			unsigned long nodesExpandedBeforeIteration = this->nodesExpanded;
			double currentDirectionBound;
			if (this->forwardSearch)
			{
				currentOriginState = originState;
				currentGoalState = goalState;
				currentDirectionBound = this->forwardBound;
			}
			else
			{
				currentOriginState = goalState;
				currentGoalState = originState;
				currentDirectionBound = this->backwardBound;
			}

			bool solved = DoIteration(currentOriginState, currentGoalState, currentOriginState, currentOriginState, currentDirectionBound, 0, midState);
			this->iteration_num++;
			unsigned long nodesExpandedThisIter = this->nodesExpanded - nodesExpandedBeforeIteration;
			if (this->verbose)
			{
				printf("Nodes expanded: %d(%d)\n", nodesExpandedThisIter, this->nodesExpanded);
			}
			if (this->forwardSearch)
			{
				this->forwardExpandedInLastIter = std::max(this->forwardExpandedInLastIter, nodesExpandedThisIter);
			}
			else
			{
				this->backwardExpandedInLastIter = std::max(this->backwardExpandedInLastIter, nodesExpandedThisIter);
			}

			if (solved)
			{
				this->pathLength = this->fBound;
				this->necessaryExpansions = this->nodesExpanded - this->lastIterBoundExpansions;
				return true;
			}
			if (this->listReady)
			{ // no solution found
				break;
			}
			if (!this->outOfSpace) // Memory is not exhausted
			{
				if (this->middleStates.size() == 0)
				{
					break;
				}
				else
				{
					this->listReady = true;
					delete this->previousBloomfilter;
					this->previousBloomfilter = nullptr;
				}
			}
			else
			{
				double saturation = this->currentBloomfilter->GetSaturation();
				if (!this->firstRun && saturation >= last_saturation)
				{
					saturationIncreasedCount += 1;
				}
				if (saturation == 1 || saturationIncreasedCount >= SATURATION_MAX_INC)
				{
					return false; //bloomfilter is fluded.
				}
				last_saturation = saturation;
				delete this->previousBloomfilter;
				this->previousBloomfilter = nullptr;
				this->listReady = false;
			}
			this->firstRun = false;
			this->forwardSearch = !this->forwardSearch;
			if (this->isConsistent)
			{
				if (this->minCurrentError != std::numeric_limits<double>::max())
				{
					this->minPreviousError = this->minCurrentError;
				}
				else
				{
					this->minPreviousError = 0;
				}
				this->minCurrentError = std::numeric_limits<double>::max();
			}
		}
		updateBounds();
	}
	return false;
}

template <class state, class action, class environment>
bool BFBDS<state, action, environment>::DoIteration(state &originState, state &goalState, state &parent, state &currState, double currentDirectionBound, double g, state &midState)
{
	this->nodesTouched++;
	double h = this->env->HCost(currState, goalState);
	if (this->isConsistent)
	{
		h += this->minPreviousError;
	}

	if (fgreater(g + h, this->fBound))
	{
		UpdateNextBound(g + h);
		return false;
	}

	if ((!this->firstRun && g + this->previousMaxG >= fBound) || g >= currentDirectionBound)
	{
		if (checkState(currState, g))
		{
			midState = currState;
			return true;
		}
		else if (g >= currentDirectionBound)
		{
			if (this->isConsistent)
			{
				this->minCurrentError = std::min(this->minCurrentError, g - this->env->HCost(currState, originState));
			}
			return false;
		}
	}

	std::map<state, double, StatesComparator> ignoreList;
	if (this->revAlgo && this->listReady)
	{
		ignoreList = std::map<state, double, StatesComparator>(this->statesComparator);
		for (auto it = this->middleStates.cbegin(); it != this->middleStates.cend();)
		{
			state perimeterState = it->first;
			double perimeterStateG = it->second;
			double calculatedF = g + this->env->HCost(currState, perimeterState) + perimeterStateG;
			if (fgreater(calculatedF, this->fBound))
			{
				UpdateNextBound(calculatedF);
				ignoreList[perimeterState] = perimeterStateG;
				it = this->middleStates.erase(it);
			}
			else
			{
				++it;
			}
		}
		if (this->middleStates.empty())
		{
			this->middleStates.insert(ignoreList.begin(), ignoreList.end());
			ignoreList.clear();
			return false;
		}
	}

	std::vector<state> neighbors;
	this->env->GetSuccessors(currState, neighbors);
	this->nodesExpanded++;
	if (g + h == this->fBound)
	{
		this->lastIterBoundExpansions++;
	}
	for (unsigned int x = 0; x < neighbors.size(); x++)
	{
		if (neighbors[x] == parent)
		{
			continue;
		}
		double edgeCost = this->env->GCost(currState, neighbors[x]);
		if (DoIteration(originState, goalState, currState, neighbors[x], currentDirectionBound, g + edgeCost, midState))
		{
			return true;
		}
	}

	if (this->revAlgo && this->listReady)
	{
		this->middleStates.insert(ignoreList.begin(), ignoreList.end());
		ignoreList.clear();
	}

	return false;
}

template <class state, class action, class environment>
bool BFBDS<state, action, environment>::checkState(state &midState, double g)
{
	if (this->listReady)
	{
		auto it = this->middleStates.find(midState);
		if (it != this->middleStates.end())
		{
			if (g + it->second <= this->fBound)
			{
				return true;
			}
		}
	}
	else if (this->firstRun || this->previousBloomfilter->Contains(this->env->GetStateHash(midState)))
	{
		this->currentMaxG = std::max(this->currentMaxG, g);
		if (this->outOfSpace)
		{
			this->currentBloomfilter->Insert(this->env->GetStateHash(midState));
		}
		else
		{
			if (this->middleStates.size() >= (int)(this->statesQuantityBound / 2))
			{
				this->outOfSpace = true;
				this->currentBloomfilter = getNextBloomfilter();
				this->currentBloomfilter->Insert(this->env->GetStateHash(midState));
				for (const auto &possibleMidStatePair : this->middleStates)
				{
					this->currentBloomfilter->Insert(this->env->GetStateHash(possibleMidStatePair.first));
				}
				this->middleStates.clear();
			}
			else
			{
				this->middleStates[midState] = g;
			}
		}
	}
	return false;
}

template <class state, class action, class environment>
void BFBDS<state, action, environment>::UpdateNextBound(double potentialNextBound)
{
	if (!fgreater(this->nextBound, this->fBound))
	{
		this->nextBound = std::max(this->fBound, potentialNextBound);
	}
	else if (fgreater(potentialNextBound, this->fBound))
	{
		this->nextBound = std::min(potentialNextBound, this->nextBound);
	}
}

template <class state, class action, class environment>
double BFBDS<state, action, environment>::calculateNextForwardBound()
{
	if (this->isUpdateByWorkload)
	{
		return calculateBoundByWorkload(this->nextBound, this->fBound, this->forwardBound, this->forwardExpandedInLastIter, this->backwardExpandedInLastIter);
	}
	return this->nextBound / 2;
}

template <class state, class action, class environment>
double BFBDS<state, action, environment>::isNextForwardSearch()
{
	return this->forwardBound >= this->backwardBound;
}

template <class state, class action, class environment>
void BFBDS<state, action, environment>::updateBounds()
{
	this->forwardBound = calculateNextForwardBound();
	this->backwardBound = this->nextBound - this->forwardBound;
	this->fBound = this->nextBound;
}

template <class state, class action, class environment>
void BFBDS<state, action, environment>::resetBfSearchValues()
{
	deleteBloomFilters();
	this->firstRun = true;
	this->forwardSearch = isNextForwardSearch();
	this->listReady = false;
	this->middleStates.clear();
	this->previousMaxG = 0;
	this->forwardExpandedInLastIter = 0;
	this->backwardExpandedInLastIter = 0;
	if (this->isConsistent)
	{
		this->minCurrentError = std::numeric_limits<double>::max();
		this->minPreviousError = 0;
	}
}

template <class state, class action, class environment>
void BFBDS<state, action, environment>::deleteBloomFilters()
{
	delete this->previousBloomfilter;
	this->previousBloomfilter = nullptr;

	delete this->currentBloomfilter;
	this->currentBloomfilter = nullptr;
}

template <class state, class action, class environment>
void BFBDS<state, action, environment>::printIterationInfo()
{
	printf("Starting iteration number: %d.\n", this->iteration_num);
	printf("Bounds: %f and %f\n", this->forwardBound, this->backwardBound);
	printf("start middleStates.size() = %d\n", this->middleStates.size());
	if (this->previousBloomfilter != nullptr)
	{
		printf("start previousBloomfilter->GetSaturation() = %f\n", this->previousBloomfilter->GetSaturation());
	}
	printf("start listReady = %d\n", this->listReady);
	printf("start firstRun = %d\n", this->firstRun);
	printf("start forwardSearch = %d\n", this->forwardSearch);
}

template <class state, class action, class environment>
int BFBDS<state, action, environment>::calculateOptimalHashNum(uint64_t bloomfilterSize, uint64_t expectedEntries)
{
	return (LN2 * bloomfilterSize) / expectedEntries;
}

template <class state, class action, class environment>
BloomFilter *BFBDS<state, action, environment>::getNextBloomfilter()
{
	int nextStartingHash = 0;
	int hashNum = DEFAULT_HASH_NUM;
	uint64_t bloomfilterSize = 0.5 * this->memoryBound;
	if (this->previousBloomfilter != nullptr)
	{
		uint64_t expectedEntries = this->previousBloomfilter->GetEntries() * pow(this->previousBloomfilter->GetSaturation(), this->previousBloomfilter->GetNumHash());
		int optimalHashNum = calculateOptimalHashNum(bloomfilterSize, expectedEntries);
		hashNum = std::max(MIN_HASH_NUM, std::min(MAX_HASH_NUM, optimalHashNum));
		nextStartingHash = this->previousBloomfilter->GetNextStartingHash();
	}
	return new BloomFilter(bloomfilterSize, hashNum, false, false, nextStartingHash);
}

#endif