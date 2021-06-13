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

template <class state, class action, class environment, class BloomFilter, bool verbose = false>
class BFBDS
{
public:
	BFBDS(environment *env, unsigned long statesQuantityBound, bool isUpdateByWorkload = true, bool isConsistent = true, bool revAlgo = false, bool F2Fheuristics = true)
	{
		this->env = env;
		this->statesQuantityBound = statesQuantityBound;
		this->isUpdateByWorkload = isUpdateByWorkload;
		this->isConsistent = isConsistent;
		this->revAlgo = revAlgo;
		this->F2Fheuristics = F2Fheuristics;
	}
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
	unsigned long memoryBound;
	BloomFilter previousBloomfilter;
	BloomFilter currentBloomfilter;
	std::vector<state> middleStates;
	bool forwardSearch;
	bool listReady;
	bool outOfSpace;
	bool firstRun;
	double backwardBound, forwardBound, fBound, nextBound;
	unsigned long lastIterBoundExpansions;
	unsigned long forwardExpandedInLastIter;
	unsigned long backwardExpandedInLastIter;
	int saturationMaxIncreasements = 5;
	double minCurrentError = std::numeric_limits<double>::max();
	double minPreviousError = 0;
};

template <class state, class action, class environment, class BloomFilter, bool verbose>
bool BFBDS<state, action, environment, BloomFilter, verbose>::GetMidState(state originState, state goalState, state &midState, std::vector<state> &thePath, int secondsLimit, bool threePhase)
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

template <class state, class action, class environment, class BloomFilter, bool verbose>
double BFBDS<state, action, environment, BloomFilter, verbose>::calculateBoundByWorkload(double newBound, double previousBound, double currentDirectionBound, uint64_t currentLoad, uint64_t currentOppositeLoad)
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

template <class state, class action, class environment, class BloomFilter, bool verbose>
bool BFBDS<state, action, environment, BloomFilter, verbose>::GetMidState(state originState, state goalState, state &midState, int secondsLimit, double startingFBound)
{
	state currentOriginState, currentGoalState;
	this->memoryBound = sizeof(originState) * statesQuantityBound;
	if (verbose)
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
			auto currentTime = std::chrono::steady_clock::now();
			std::chrono::duration<double> elapsed_seconds = currentTime - startTime;
			if (elapsed_seconds.count() >= secondsLimit)
			{
				return false;
			}
			if (verbose)
			{
				printf("Bounds: %f and %f\n", this->forwardBound, this->backwardBound);
				printf("Starting iteration number: %d.\n", this->iteration_num);
				printf("start middleStates.size() = %d\n", this->middleStates.size());
				printf("start previousBloomfilter.getSaturationgetSaturation() = %f\n", this->previousBloomfilter.getSaturation());
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
			if (verbose)
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
			if (!this->outOfSpace)
			{
				if (this->middleStates.size() == 0)
				{
					break;
				}
				else
				{
					this->listReady = true;
					this->previousBloomfilter.clear();
				}
			}
			else
			{
				double saturation = this->currentBloomfilter.getSaturation();
				if (!this->firstRun && saturation >= last_saturation)
				{
					saturationIncreasedCount += 1;
				}
				if (saturation == 1 || saturationIncreasedCount >= this->saturationMaxIncreasements)
				{
					return false; //bloomfilter is fluded.
				}
				last_saturation = saturation;
				this->previousBloomfilter.clear();
				this->previousBloomfilter = currentBloomfilter;
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
	return false;
}

template <class state, class action, class environment, class BloomFilter, bool verbose>
bool BFBDS<state, action, environment, BloomFilter, verbose>::DoIteration(state &originState, state &goalState, state &parent, state &currState, double currentDirectionBound, double g, state &midState)
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

template <class state, class action, class environment, class BloomFilter, bool verbose>
bool BFBDS<state, action, environment, BloomFilter, verbose>::checkState(state &midState)
{
	if (this->listReady)
	{
		if (std::find(this->middleStates.begin(), this->middleStates.end(), midState) != this->middleStates.end())
		{
			return true;
		}
	}
	else if (this->firstRun || this->previousBloomfilter.contains(midState))
	{
		if (this->outOfSpace)
		{
			this->currentBloomfilter.insert(midState);
		}
		else
		{
			if (this->middleStates.size() >= (int)(this->statesQuantityBound / 2))
			{
				this->outOfSpace = true;
				if (this->previousBloomfilter.getCount() > 0)
				{
					unsigned int nextK = std::max(1, std::min(10, int(0.693 * (this->memoryBound / 2) / (this->previousBloomfilter.getCount() * this->previousBloomfilter.getSaturation()))));
					this->currentBloomfilter = BloomFilter(this->memoryBound / 2, nextK, this->previousBloomfilter.getNextOffset());
				}
				else
				{
					this->currentBloomfilter = BloomFilter(this->memoryBound / 2, 5, this->previousBloomfilter.getNextOffset());
				}
				for (state possibleMidState : this->middleStates)
				{
					this->currentBloomfilter.insert(possibleMidState);
				}
				this->currentBloomfilter.insert(midState);
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

template <class state, class action, class environment, class BloomFilter, bool verbose>
void BFBDS<state, action, environment, BloomFilter, verbose>::UpdateNextBound(double potentialNextBound)
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

template <class state, class action, class environment, class BloomFilter, bool verbose>
double BFBDS<state, action, environment, BloomFilter, verbose>::calculateNextForwardBound()
{
	if (this->isUpdateByWorkload)
	{
		return calculateBoundByWorkload(this->nextBound, this->fBound, this->forwardBound, this->forwardExpandedInLastIter, this->backwardExpandedInLastIter);
	}
	return ceil(this->nextBound / 2);
}

template <class state, class action, class environment, class BloomFilter, bool verbose>
void BFBDS<state, action, environment, BloomFilter, verbose>::updateBounds()
{
	this->forwardBound = calculateNextForwardBound();
	this->fBound = this->nextBound;
	this->backwardBound = this->fBound - this->forwardBound;
}

template <class state, class action, class environment, class BloomFilter, bool verbose>
void BFBDS<state, action, environment, BloomFilter, verbose>::resetBfSearchValues()
{
	this->previousBloomfilter.clear();
	this->currentBloomfilter.clear();
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

#endif