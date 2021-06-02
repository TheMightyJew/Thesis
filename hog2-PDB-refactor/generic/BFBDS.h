/*
 *  BFBDS.h
 *
 *  Created by Steven Danishevski on DEC 19.
 *
 */

#ifndef BFBDS_H
#define BFBDS_H

#include <iostream>
#include <algorithm>
#include <unordered_set>
#include "SearchEnvironment.h"
#include <math.h>
#include "MM.h"
#include "IDBiHS.h"

template <class state, class action, class environment, class BloomFilter, bool verbose = false>
class BFBDS
{
public:
	BFBDS(unsigned long statesQuantityBound, bool isUpdateByWorkload = true, bool isConsistent = true, bool revAlgo = false, bool F2Fheuristics = true)
	{
		this->statesQuantityBound = statesQuantityBound;
		this->isUpdateByWorkload = isUpdateByWorkload;
		this->isConsistent = isConsistent;
		this->revAlgo = revAlgo;
		this->F2Fheuristics = F2Fheuristics;
	}
	virtual ~BFBDS() {}
	bool GetMidState(environment *env, state fromState, state toState, state &midState, std::vector<state> &thePath, int secondsLimit = 600, bool threePhase = true);
	double getPathLength() { return pathLength; }
	uint64_t getNodesExpanded() { return nodesExpanded; }
	uint64_t getNecessaryExpansions() { return necessaryExpansions; }
	int getIterationsNum() { return iteration_num; }
	void ResetNodeCount() { nodesExpanded = nodesTouched = 0; }
	bool isThreePhase() { return threePhase; }

private:
	uint64_t nodesExpanded, nodesTouched, statesQuantityBound, necessaryExpansions;
	int iteration_num;
	double pathLength;
	unsigned long memoryBound;
	bool threePhase;
	bool isUpdateByWorkload;
	bool isConsistent;
	bool revAlgo;
	bool F2Fheuristics;

	bool GetMidState(environment *env, state fromState, state toState, state &midState, int secondsLimit = 600, double startingFBound = 0);
	void UpdateNextBound(double fCost);
	double updateBoundByWorkload(double newbound, double prevBound, double oldForwardBound, uint64_t forwardLoad, uint64_t backwardLoad);
	bool DoIteration(SearchEnvironment<state, action> *env, state &parent, state &currState, double bound, double g, state &midState);
	bool checkState(state &midState);

	BloomFilter previousBloomfilter;
	BloomFilter currentBloomfilter;
	bool forwardSearch;
	bool listReady;
	bool outOfSpace;
	bool firstRun;
	double backwardBound, forwardBound, fBound, nextBound;
	unsigned long lastIterBoundExpansions;
	unsigned long forwardExpandedInLastIter = 0;
	unsigned long backwardExpandedInLastIter = 0;
	unsigned long lastBLinsertions = 0;
	std::vector<state> middleStates;
	double minCurrentError = std::numeric_limits<double>::max();
	double minPreviousError = 0;
	state goal;
	state from;
};

template <class state, class action, class environment, class BloomFilter, bool verbose>
bool BFBDS<state, action, environment, BloomFilter, verbose>::GetMidState(environment *env,
																		  state fromState, state toState, state &midState, std::vector<state> &thePath, int secondsLimit, bool threePhase)
{
	bool solved = false;
	nodesExpanded = 0;
	nodesTouched = 0;
	necessaryExpansions = 0;
	unsigned long currentNodesExapanded = 0;
	double lastBound = 0;
	if (threePhase)
	{
		MM<state, action, environment> mm;
		solved = mm.GetPath(env, fromState, toState, env, env, thePath, secondsLimit, statesQuantityBound);
		currentNodesExapanded += mm.GetNodesExpanded();
		necessaryExpansions = mm.GetNecessaryExpansions();
		lastBound = mm.getLastBound();
	}
	if (solved)
	{
		pathLength = env->GetPathLength(thePath);
		nodesExpanded = currentNodesExapanded;
		return true;
	}
	else
	{
		solved = GetMidState(env, fromState, toState, midState, secondsLimit, lastBound);
		currentNodesExapanded += nodesExpanded;
		lastBound = backwardBound + forwardBound;
		if (solved)
		{
			nodesExpanded = currentNodesExapanded;
			return true;
		}
		else
		{
			IDBiHS<environment, state, action, false> idbihs(env, F2Fheuristics, isConsistent, isUpdateByWorkload);
			solved = idbihs.GetMidState(fromState, toState, midState, secondsLimit, lastBound);
			currentNodesExapanded += idbihs.GetNodesExpanded();
			if (solved)
			{
				pathLength = idbihs.getPathLength();
				nodesExpanded = currentNodesExapanded;
				necessaryExpansions = idbihs.GetNecessaryExpansions();
				return true;
			}
			else
			{
				return false;
			}
		}
	}
}

template <class state, class action, class environment, class BloomFilter, bool verbose>
double BFBDS<state, action, environment, BloomFilter, verbose>::updateBoundByWorkload(double newbound, double prevBound, double oldForwardBound, uint64_t forwardLoad, uint64_t backwardLoad)
{
	if (forwardLoad <= backwardLoad)
	{
		return oldForwardBound + newbound - prevBound;
	}
	else
	{
		return oldForwardBound;
	}
}

template <class state, class action, class environment, class BloomFilter, bool verbose>
bool BFBDS<state, action, environment, BloomFilter, verbose>::GetMidState(environment *env,
																		  state fromState, state toState, state &midState, int secondsLimit, double startingFBound)
{
	memoryBound = sizeof(fromState) * statesQuantityBound;
	if (verbose)
	{
		printf("\nmemory is %1.1llu\n", memoryBound);
		printf("\nStarting to solve with BFBDS\n");
	}
	nodesExpanded = nodesTouched = 0;
	double initialHeuristic = env->HCost(fromState, toState);
	fBound = round(std::max(initialHeuristic, startingFBound));
	nextBound = fBound;
	forwardBound = ceil(fBound / 2);
	backwardBound = fBound - forwardBound;
	int saturationIncreased = 0;
	int saturationMaxIncreasements = 10;
	iteration_num = 0;
	double bound;
	firstRun = true;
	forwardSearch = true;
	listReady = false;
	middleStates.clear();
	unsigned long nodesExpandedSoFar = 0;
	double last_saturation = 1;
	auto startTime = std::chrono::steady_clock::now();
	while (true)
	{
		lastIterBoundExpansions = 0;
		if (verbose)
		{
			printf("Bounds: %f and %f\n", forwardBound, backwardBound);
		}
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
				printf("Starting iteration number: %d. ", iteration_num);
				printf("start middleStates.size() = %d\n", middleStates.size());
				printf("start previousBloomfilter.getSaturationgetSaturation() = %f\n", previousBloomfilter.getSaturation());
				printf("start listReady = %d\n", listReady);
				printf("start firstRun = %d\n", firstRun);
				printf("start forwardSearch = %d\n", forwardSearch);
			}
			if (forwardSearch)
			{
				forwardExpandedInLastIter = 0;
				bound = forwardBound;
				from = fromState;
				goal = toState;
			}
			else
			{
				backwardExpandedInLastIter = 0;
				bound = backwardBound;
				from = toState;
				goal = fromState;
			}
			outOfSpace = false;
			if (revAlgo && listReady)
			{
				for (uint64_t i = middleStates.size() - 1; i != static_cast<uint64_t>(-1); i--)
				{
					state perimeterState = middleStates[i];
					double calculatedF = env->HCost(from, perimeterState) + (fBound - bound);
					if (fgreater(calculatedF, fBound))
					{
						UpdateNextBound(calculatedF);
						middleStates[i] = middleStates.back();
						middleStates.pop_back();
					}
				}
			}
			bool solved = DoIteration(env, from, from, bound, 0, midState);
			unsigned long nodesExpandedThisIter = nodesExpanded - nodesExpandedSoFar;
			nodesExpandedSoFar = nodesExpanded;
			if (verbose)
			{
				printf("Nodes expanded: %d(%d)\n", nodesExpandedThisIter, nodesExpanded);
			}
			if (forwardSearch)
			{
				forwardExpandedInLastIter = std::max(forwardExpandedInLastIter, nodesExpandedThisIter);
			}
			else
			{
				backwardExpandedInLastIter = std::max(backwardExpandedInLastIter, nodesExpandedThisIter);
				;
			}
			iteration_num++;
			if (solved)
			{
				pathLength = fBound;
				necessaryExpansions = nodesExpanded - lastIterBoundExpansions;
				return true;
			}
			if (listReady)
			{ // no solution found
				if (verbose)
				{
					std::cout << "No solution" << std::endl;
				}
				break;
			}
			if (!outOfSpace)
			{
				if (middleStates.size() == 0)
				{
					if (verbose)
					{
						std::cout << "No solution" << std::endl;
					}
					break;
				}
				else
				{
					listReady = true;
					previousBloomfilter.clear();
				}
			}
			else
			{
				double saturation = currentBloomfilter.getSaturation();
				if (saturation >= last_saturation)
				{
					saturationIncreased += 1;
				}
				if (saturation == 1 || saturationIncreased >= saturationMaxIncreasements)
				{
					if (verbose)
					{
						std::cout << "\t\tBloomFilter Overflow" << std::endl;
					}
					return false; //bloomfilter is fluded.
				}
				last_saturation = saturation;
				if (verbose)
				{
					printf("Bloomfilter saturation is: %1.3f%%\n", saturation);
				}
				previousBloomfilter.clear();
				previousBloomfilter = currentBloomfilter;
				listReady = false;
			}
			firstRun = false;
			forwardSearch = !forwardSearch;
			if (minCurrentError != std::numeric_limits<double>::max())
			{
				minPreviousError = minCurrentError;
			}
			else
			{
				minPreviousError = 0;
			}
			minCurrentError = std::numeric_limits<double>::max();
		}

		last_saturation = 1;
		lastBLinsertions = 0;
		saturationIncreased = 0;
		previousBloomfilter.clear();
		firstRun = true;
		forwardSearch = true;
		listReady = false;
		middleStates.clear();
		minCurrentError = std::numeric_limits<double>::max();
		minPreviousError = 0;

		if (!isUpdateByWorkload)
		{
			forwardBound = ceil(nextBound / 2);
		}
		else
		{
			forwardBound = updateBoundByWorkload(nextBound, fBound, forwardBound, forwardExpandedInLastIter, backwardExpandedInLastIter);
		}
		fBound = nextBound;
		backwardBound = fBound - forwardBound;
		forwardExpandedInLastIter = 0;
		backwardExpandedInLastIter = 0;
	}
	return false;
}

template <class state, class action, class environment, class BloomFilter, bool verbose>
bool BFBDS<state, action, environment, BloomFilter, verbose>::DoIteration(SearchEnvironment<state, action> *env,
																		  state &parent, state &currState, double bound, double g, state &midState)
{
	double h = env->HCost(currState, goal);
	std::vector<state> ignoreList;
	if (isConsistent && g < bound)
	{
		h += minPreviousError;
	}

	if (g > bound)
	{
		return false;
	}

	if (fgreater(g + h, fBound))
	{
		UpdateNextBound(g + h);
		return false;
	}

	if (g == bound)
	{
		if (checkState(currState))
		{
			midState = currState;
			return true;
		}
		else
		{
			if (isConsistent)
			{
				minCurrentError = std::min(minCurrentError, g - env->HCost(currState, from));
			}
			return false;
		}
	}

	if (revAlgo && listReady)
	{
		for (uint64_t i = middleStates.size() - 1; i != static_cast<uint64_t>(-1); i--)
		{
			state perimeterState = middleStates[i];
			double calculatedF = g + env->HCost(currState, perimeterState) + (fBound - bound);
			if (fgreater(calculatedF, fBound))
			{
				UpdateNextBound(calculatedF);
				ignoreList.push_back(perimeterState);
				middleStates[i] = middleStates.back();
				middleStates.pop_back();
			}
		}
		if (middleStates.size() == 0)
		{
			middleStates.insert(middleStates.end(), ignoreList.begin(), ignoreList.end());
			ignoreList.clear();
			return false;
		}
	}

	std::vector<state> neighbors;
	env->GetSuccessors(currState, neighbors);
	nodesTouched += neighbors.size();
	nodesExpanded++;
	if (g + h == fBound)
	{
		lastIterBoundExpansions++;
	}
	for (unsigned int x = 0; x < neighbors.size(); x++)
	{
		if (neighbors[x] == parent)
		{
			continue;
		}
		double edgeCost = env->GCost(currState, neighbors[x]);
		if (DoIteration(env, currState, neighbors[x], bound, g + edgeCost, midState))
		{
			return true;
		}
	}

	if (revAlgo && listReady)
	{
		middleStates.insert(middleStates.end(), ignoreList.begin(), ignoreList.end());
		ignoreList.clear();
	}

	return false;
}

template <class state, class action, class environment, class BloomFilter, bool verbose>
bool BFBDS<state, action, environment, BloomFilter, verbose>::checkState(state &midState)
{
	if (listReady)
	{
		if (std::find(middleStates.begin(), middleStates.end(), midState) != middleStates.end())
		{
			return true;
		}
	}
	else if (firstRun || previousBloomfilter.contains(midState))
	{
		if (outOfSpace)
		{
			currentBloomfilter.insert(midState);
			lastBLinsertions += 1;
		}
		else
		{
			if (middleStates.size() >= (int)(statesQuantityBound / 2))
			{
				outOfSpace = true;
				if (lastBLinsertions > 0)
					currentBloomfilter = BloomFilter(memoryBound / 2, round(0.693 * (memoryBound / 2) / (lastBLinsertions)), previousBloomfilter.hashOffset + previousBloomfilter.getK());
				else
					currentBloomfilter = BloomFilter(memoryBound / 2, 5, previousBloomfilter.hashOffset + previousBloomfilter.getK());
				lastBLinsertions = middleStates.size() + 1;
				for (state possibleMidState : middleStates)
				{
					currentBloomfilter.insert(possibleMidState);
				}
				currentBloomfilter.insert(midState);
				middleStates.clear();
			}
			else
			{
				middleStates.push_back(midState);
			}
		}
	}
	return false;
}

template <class state, class action, class environment, class BloomFilter, bool verbose>
void BFBDS<state, action, environment, BloomFilter, verbose>::UpdateNextBound(double fCost)
{
	if (!fgreater(nextBound, fBound))
	{
		nextBound = std::max(fBound, fCost);
	}
	else if (fgreater(fCost, fBound))
	{
		nextBound = std::min(fCost, nextBound);
	}
}

#endif