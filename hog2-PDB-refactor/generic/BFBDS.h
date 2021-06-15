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

template <class state, class action, class environment>
class BFBDS
{
public:
	BFBDS(environment *env, unsigned long statesQuantityBound, bool isUpdateByWorkload = true, bool isConsistent = true, bool revAlgo = false, bool F2Fheuristics = true, bool verbose = false);
	virtual ~BFBDS() {}
	bool GetMidState(state originState, state goalState, state &midState, std::vector<state> &thePath, int secondsLimit = 600, bool threePhase = true);
	double getPathLength() { return pathLength; }
	uint64_t getNodesExpanded() { return nodesExpanded; }
	uint64_t getNecessaryExpansions() { return necessaryExpansions; }
	int getIterationsNum() { return iteration_num; }

private:
	bool GetMidState(state originState, state goalState, state &midState, int secondsLimit = 600, double startingFBound = 0);
	bool DoIteration(state &originState, state &goalState, state &parent, state &currState, double currentDirectionBound, double g, state &midState);
	bool checkState(state &midState);
	void updateBounds();
	double calculateNextForwardBound();
	void UpdateNextBound(double potentialNextBound);
	double calculateBoundByWorkload(double newBound, double previousBound, double currentDirectionBound, uint64_t currentLoad, uint64_t currentOppositeLoad);
	int calculateOptimalHashNum(uint64_t bloomfilterSize, uint64_t expectedEntries);
	void resetBfSearchValues();

	//Constructor variables
	environment *env;
	bool threePhase;
	bool isUpdateByWorkload;
	bool isConsistent;
	bool revAlgo;
	bool F2Fheuristics;
	uint64_t statesQuantityBound;

	//Results variables
	uint64_t nodesExpanded, nodesTouched, necessaryExpansions;
	int iteration_num;
	double pathLength;

	//Run variables
	bool verbose;
	unsigned long memoryBound;
	BloomFilter *previousBloomfilter = nullptr;
	BloomFilter *currentBloomfilter = nullptr;
	std::vector<state> middleStates;
	bool forwardSearch;
	bool listReady;
	bool outOfSpace;
	bool firstRun;
	double backwardBound, forwardBound, fBound, nextBound;
	unsigned long lastIterBoundExpansions;
	unsigned long forwardExpandedInLastIter;
	unsigned long backwardExpandedInLastIter;
	int SATURATION_MAX_INC = 5;
	int DEFAULT_HASH_NUM = 5;
	int MIN_HASH_NUM = 1;
	int MAX_HASH_NUM = 10;
	double LN2 = 0.693;
	double minCurrentError = std::numeric_limits<double>::max();
	double minPreviousError = 0;
};

template <class state, class action, class environment>
BFBDS<state, action, environment>::BFBDS(environment *env, unsigned long statesQuantityBound, bool isUpdateByWorkload, bool isConsistent, bool revAlgo, bool F2Fheuristics, bool verbose) : env(env), statesQuantityBound(statesQuantityBound), isUpdateByWorkload(isUpdateByWorkload), isConsistent(isConsistent), revAlgo(revAlgo), F2Fheuristics(F2Fheuristics), verbose(verbose) {}

template <class state, class action, class environment>
bool BFBDS<state, action, environment>::GetMidState(state originState, state goalState, state &midState, std::vector<state> &thePath, int secondsLimit, bool threePhase)
{
	this->nodesExpanded = this->nodesTouched = this->necessaryExpansions = 0;
	double lastBound;
	auto startTime = std::chrono::steady_clock::now();
	double elapsed_seconds;
	if (threePhase)
	{
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
	elapsed_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();
	bool solved = GetMidState(originState, goalState, midState, secondsLimit - elapsed_seconds, lastBound);
	lastBound = this->backwardBound + this->forwardBound;
	if (solved)
	{
		return true;
	}
	else
	{
		IDBiHS<environment, state, action, false> idbihs(this->env, this->F2Fheuristics, this->isConsistent, this->isUpdateByWorkload);
		elapsed_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();
		solved = idbihs.GetMidState(originState, goalState, midState, secondsLimit - elapsed_seconds, lastBound);
		unsigned long test = this->nodesExpanded;
		this->nodesExpanded += idbihs.GetNodesExpanded();
		if (solved)
		{
			this->pathLength = idbihs.getPathLength();
			this->necessaryExpansions = idbihs.GetNecessaryExpansions();
		}
		return solved;
	}
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
	if (this->verbose)
	{
		printf("\nmemory is %1.1llu\n", this->memoryBound);
		printf("\nStarting to solve with BFBDS\n");
	}
	double initialHeuristic = this->env->HCost(originState, goalState);
	this->nextBound = this->fBound = round(std::max(initialHeuristic, startingFBound));
	this->forwardBound = ceil(this->fBound / 2);
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
			this->currentBloomfilter = nullptr;
			auto currentTime = std::chrono::steady_clock::now();
			std::chrono::duration<double> elapsed_seconds = currentTime - startTime;
			if (elapsed_seconds.count() >= secondsLimit)
			{
				resetBfSearchValues();
				return false;
			}
			if (this->verbose)
			{
				printf("Bounds: %f and %f\n", this->forwardBound, this->backwardBound);
				printf("Starting iteration number: %d.\n", this->iteration_num);
				printf("start middleStates.size() = %d\n", this->middleStates.size());
				if (this->previousBloomfilter != nullptr)
				{
					printf("start previousBloomfilter->GetSaturation() = %f\n", this->previousBloomfilter->GetSaturation());
				}
				printf("start listReady = %d\n", this->listReady);
				printf("start firstRun = %d\n", this->firstRun);
				printf("start forwardSearch = %d\n", this->forwardSearch);
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
				resetBfSearchValues();
				return true;
			}
			if (this->listReady)
			{ // no solution found
				break;
			}
			if (!this->outOfSpace)
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
				if (saturation == 1 || saturationIncreasedCount >= this->SATURATION_MAX_INC)
				{
					resetBfSearchValues();
					return false; //bloomfilter is fluded.
				}
				last_saturation = saturation;
				delete this->previousBloomfilter;
				this->previousBloomfilter = nullptr;
				this->previousBloomfilter = this->currentBloomfilter;
				this->listReady = false;
			}
			this->firstRun = false;
			this->forwardSearch = !this->forwardSearch;
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
		updateBounds();
	}
	resetBfSearchValues();
	return false;
}

template <class state, class action, class environment>
bool BFBDS<state, action, environment>::DoIteration(state &originState, state &goalState, state &parent, state &currState, double currentDirectionBound, double g, state &midState)
{
	double h = this->env->HCost(currState, goalState);
	std::vector<state> ignoreList;
	if (this->isConsistent && g < currentDirectionBound)
	{
		h += this->minPreviousError;
	}

	if (g > currentDirectionBound)
	{
		return false;
	}

	if (fgreater(g + h, this->fBound))
	{
		UpdateNextBound(g + h);
		return false;
	}

	if (g == currentDirectionBound)
	{
		if (checkState(currState))
		{
			midState = currState;
			return true;
		}
		else
		{
			if (this->isConsistent)
			{
				this->minCurrentError = std::min(this->minCurrentError, g - this->env->HCost(currState, originState));
			}
			return false;
		}
	}

	if (this->revAlgo && this->listReady)
	{
		for (auto it = this->middleStates.rbegin(); it != this->middleStates.rend(); ++it)
		{
			state perimeterState = *it;
			double calculatedF = g + this->env->HCost(currState, perimeterState) + (this->fBound - currentDirectionBound);
			if (fgreater(calculatedF, this->fBound))
			{
				UpdateNextBound(calculatedF);
				ignoreList.push_back(perimeterState);
				this->middleStates.erase(std::next(it).base()); //The correct way to erase with rev-iterator
			}
		}
		if (this->middleStates.empty())
		{
			this->middleStates.insert(this->middleStates.end(), ignoreList.begin(), ignoreList.end());
			ignoreList.clear();
			return false;
		}
	}

	std::vector<state> neighbors;
	this->env->GetSuccessors(currState, neighbors);
	this->nodesTouched += neighbors.size();
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
		this->middleStates.insert(this->middleStates.end(), ignoreList.begin(), ignoreList.end());
		ignoreList.clear();
	}

	return false;
}

template <class state, class action, class environment>
bool BFBDS<state, action, environment>::checkState(state &midState)
{
	if (this->listReady)
	{
		if (std::find(this->middleStates.begin(), this->middleStates.end(), midState) != this->middleStates.end())
		{
			return true;
		}
	}
	else if (this->firstRun || this->previousBloomfilter->Contains(this->env->GetStateHash(midState)))
	{
		if (this->outOfSpace)
		{
			this->currentBloomfilter->Insert(this->env->GetStateHash(midState));
		}
		else
		{
			if (this->middleStates.size() > (int)(this->statesQuantityBound / 2))
			{
				this->outOfSpace = true;
				int nextStartingHash = 0;
				int hashNum = DEFAULT_HASH_NUM;
				uint64_t bloomfilterSize = 0.5 * this->memoryBound;
				if (this->previousBloomfilter != nullptr)
				{
					uint64_t expectedEntries = this->previousBloomfilter->GetEntries() * this->previousBloomfilter->GetSaturation();
					int optimalHashNum = calculateOptimalHashNum(bloomfilterSize, expectedEntries);
					hashNum = std::max(MIN_HASH_NUM, std::min(MAX_HASH_NUM, optimalHashNum));
					nextStartingHash = this->previousBloomfilter->GetNextStartingHash();
				}
				this->currentBloomfilter = new BloomFilter(static_cast<uint64_t>(bloomfilterSize), hashNum, false, false, nextStartingHash);

				this->currentBloomfilter->Insert(this->env->GetStateHash(midState));
				for (state possibleMidState : this->middleStates)
				{
					this->currentBloomfilter->Insert(this->env->GetStateHash(possibleMidState));
				}
				
				this->middleStates.clear();
			}
			else
			{
				this->middleStates.push_back(midState);
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
		this->nextBound = std::max(fBound, potentialNextBound);
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
	return ceil(this->nextBound / 2);
}

template <class state, class action, class environment>
void BFBDS<state, action, environment>::updateBounds()
{
	this->forwardBound = calculateNextForwardBound();
	this->fBound = this->nextBound;
	this->backwardBound = this->fBound - this->forwardBound;
}

template <class state, class action, class environment>
void BFBDS<state, action, environment>::resetBfSearchValues()
{
	delete this->previousBloomfilter;
	this->previousBloomfilter = nullptr;
	delete this->currentBloomfilter;
	this->currentBloomfilter = nullptr;
	this->middleStates.clear();
	this->firstRun = true;
	this->forwardSearch = true;
	this->listReady = false;
	this->middleStates.clear();
	this->minCurrentError = std::numeric_limits<double>::max();
	this->minPreviousError = 0;
	this->forwardExpandedInLastIter = 0;
	this->backwardExpandedInLastIter = 0;
}

template <class state, class action, class environment>
int BFBDS<state, action, environment>::calculateOptimalHashNum(uint64_t bloomfilterSize, uint64_t expectedEntries){
	return (this->LN2 * bloomfilterSize) / expectedEntries;
}

#endif