/*
 *  MBBDS.h
 *
 *  Created by Steven Danishevski on DEC 19.
 *
 */

#ifndef MBBDS_H
#define MBBDS_H

#include <iostream>
#include<algorithm> 
#include <unordered_set>
#include "SearchEnvironment.h"
#include <math.h>

template <class state, class action, class BloomFilter, bool verbose = true>
class MBBDS {
public:
	MBBDS(unsigned long statesQuantityBound, bool isUpdateByWorkload=false, bool isConsistent = false) {
		this->statesQuantityBound = statesQuantityBound;
		this->isUpdateByWorkload = isUpdateByWorkload;
		this->isConsistent = isConsistent;
	}
	virtual ~MBBDS() {}
	bool GetMidState(SearchEnvironment<state, action>* env, state fromState, state toState, state &midState, int secondsLimit=600, double startingFBound=0);
	double getPathLength()	{ return pathLength; }
	uint64_t GetNodesExpanded() { return nodesExpanded; }
	uint64_t GetNodesTouched() { return nodesTouched; }
	uint64_t GetNecessaryExpansions() { return necessaryExpansions; }
	double getLastBound() { return backwardBound+forwardBound; }
	int getIterationNum() { return iteration_num; }
	void ResetNodeCount() { nodesExpanded = nodesTouched = 0; }
private:
	void UpdateNextBound(double fCost);
	double updateBoundByWorkload(double newbound, double prevBound, double oldForwardBound,uint64_t forwardLoad,uint64_t backwardLoad);
	unsigned long nodesExpanded, nodesTouched, statesQuantityBound, necessaryExpansions, lastIterBoundExpansions;
;
	int iteration_num;
	BloomFilter previousBloomfilter;
	BloomFilter currentBloomfilter;
	bool listReady;
	bool outOfSpace;
	bool firstRun;
	double backwardBound, forwardBound , fBound;
	unsigned long forwardExpandedInLastIter = 0;
	unsigned long backwardExpandedInLastIter = 0;
	double pathLength;
	std::unordered_set<state> middleStates;
	unsigned long memoryBound;

	bool DoIteration(SearchEnvironment<state, action>* env,
		state parent, state currState, double bound, double g, state& midState);
	bool checkState(state midState);
	double nextBound;
	bool forwardSearch;
	bool isUpdateByWorkload;
	bool isConsistent;
	double minCurrentError = std::numeric_limits<double>::max();
	double minPreviousError = 0;
	state goal;
	state from;



};
template <class state, class action, class BloomFilter, bool verbose>
double MBBDS<state, action, BloomFilter, verbose>::updateBoundByWorkload(double newbound, double prevBound, double oldForwardBound, uint64_t forwardLoad,uint64_t backwardLoad){
	//printf("forwardLoad: %llu,backwardLoad: %llu,oldForwardBound%d\n",forwardLoad,backwardLoad,oldForwardBound);
	if (forwardLoad <= backwardLoad){
		return oldForwardBound + newbound - prevBound;
	}
	else{
		return oldForwardBound;
	}
}

template <class state, class action, class BloomFilter, bool verbose>
bool MBBDS<state, action, BloomFilter, verbose>::GetMidState(SearchEnvironment<state, action>* env,
	state fromState, state toState, state &midState, int secondsLimit, double startingFBound)
{
	memoryBound = sizeof(fromState)*statesQuantityBound;
	if(verbose){
		printf("\nmemory is %1.1llu\n", memoryBound);
		printf("\nStarting to solve with MBBDS\n");
	}
	nodesExpanded = nodesTouched = 0;
	double initialHeuristic = env->HCost(fromState, toState);
	fBound = round(std::max(initialHeuristic, startingFBound));
	nextBound = fBound;
	forwardBound = ceil(startingFBound / 2);
	backwardBound = startingFBound - forwardBound;
	int saturationIncreased = 0;
	//changed
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
	while(true){
		lastIterBoundExpansions = 0;
		if (verbose){
			printf("Bounds: %f and %f\n", forwardBound, backwardBound);
		}
		while (true){
			auto currentTime = std::chrono::steady_clock::now();
			std::chrono::duration<double> elapsed_seconds = currentTime-startTime;
			if(elapsed_seconds.count() >= secondsLimit){
				return false;
			}
			if (verbose){
				printf("Starting iteration number: %d. ", iteration_num);
				printf("start middleStates.size() = %d\n", middleStates.size());
				printf("start previousBloomfilter.getSaturationgetSaturation() = %f\n", previousBloomfilter.getSaturation());
				printf("start listReady = %d\n", listReady);
				printf("start firstRun = %d\n", firstRun);
				printf("start forwardSearch = %d\n", forwardSearch);
			}
			if (forwardSearch) {
				forwardExpandedInLastIter = 0;
				bound = forwardBound;
				from = fromState;
				goal = toState;
			}
			else {
				backwardExpandedInLastIter = 0;
				bound = backwardBound;
				from = toState;
				goal = fromState;
			}
			outOfSpace = false;
			bool solved = DoIteration(env, from, from, bound, 0, midState);
			unsigned long nodesExpandedThisIter = nodesExpanded-nodesExpandedSoFar;
			nodesExpandedSoFar = nodesExpanded;
			if(verbose){
				printf("Nodes expanded: %d(%d)\n", nodesExpandedThisIter, nodesExpanded);
			}
			if (forwardSearch) {
				forwardExpandedInLastIter = std::max(forwardExpandedInLastIter, nodesExpandedThisIter);
			}
			else {
				backwardExpandedInLastIter = std::max(backwardExpandedInLastIter, nodesExpandedThisIter);;
			}
			iteration_num++;
			if (solved) {
				pathLength = fBound;
				necessaryExpansions = nodesExpanded - lastIterBoundExpansions;
				return true;
			}
			if (listReady) {// no solution found
				if(verbose){
					std::cout << "No solution" << std::endl;
				}
				break;
			}
			if (!outOfSpace) {
				if(middleStates.size()==0){
					if(verbose){
						std::cout << "No solution" << std::endl;
					}
					break;
				}
				else{
					listReady = true;
					previousBloomfilter.clear();
				}				
			}
			else {
				double saturation = currentBloomfilter.getSaturation();
				if(saturation >= last_saturation){
					saturationIncreased += 1;
				}
				if (saturation == 1 || saturationIncreased >= saturationMaxIncreasements) {
					if(verbose){
						std::cout << "\t\tBloomFilter Overflow" << std::endl;
					}
					return false; //bloomfilter is fluded.
				}
				last_saturation = saturation;
				if(verbose){
					printf("Bloomfilter saturation is: %1.3f%%\n", saturation);
				}
				previousBloomfilter.clear();
				previousBloomfilter = currentBloomfilter;
				listReady = false;
			}
			firstRun = false;
			forwardSearch = !forwardSearch;
			if(minCurrentError != std::numeric_limits<double>::max()){
				minPreviousError = minCurrentError;
			}
			else{
				minPreviousError = 0;
			}
			minCurrentError = std::numeric_limits<double>::max();
		}	
		
		last_saturation = 1;
		saturationIncreased = 0;
		previousBloomfilter.clear();
		firstRun = true;
		forwardSearch = true;
		listReady = false;
		middleStates.clear();
		minCurrentError = std::numeric_limits<double>::max();
		minPreviousError = 0;
		nextBound = std::max(nextBound, forwardBound+backwardBound+1);
		
		if (!isUpdateByWorkload){
			forwardBound = ceil(nextBound / 2);
		} 
		else{
			forwardBound = updateBoundByWorkload(nextBound, fBound, forwardBound, forwardExpandedInLastIter, backwardExpandedInLastIter);
		}
		fBound = nextBound;	
		backwardBound = fBound - forwardBound;
		forwardExpandedInLastIter = 0;
		backwardExpandedInLastIter = 0;
	}
	return false;
}


template <class state, class action, class BloomFilter, bool verbose>
bool MBBDS<state, action, BloomFilter, verbose>::DoIteration(SearchEnvironment<state, action>* env,
	state parent, state currState, double bound, double g, state& midState)
{
	double h = env->HCost(currState, goal);
	
	if(isConsistent && g < bound){
		h += minPreviousError;
	}
	
	if (g > bound  || fgreater(g + h, fBound))
	{
		UpdateNextBound(g + h);
		return false;
	}

	if (g == bound) {
    if (isConsistent){
      minCurrentError = std::min(minCurrentError, g - env->HCost(currState,from)); 
    }
		if (checkState(currState)){
			midState = currState;
			return true;
		}
		else{
			return false;
		}
	}
	
	std::vector<state> neighbors;
	env->GetSuccessors(currState, neighbors);
	nodesTouched += neighbors.size();
	nodesExpanded++;
	if(g + h == fBound){
		lastIterBoundExpansions++;
	}
	for (unsigned int x = 0; x < neighbors.size(); x++)
	{
		if (neighbors[x] == parent) {
			continue;
		}
		double edgeCost = env->GCost(currState, neighbors[x]);
		if (DoIteration(env, currState, neighbors[x], bound, g + edgeCost, midState)) {
			return true;
		}
	}
	return false;
}

template<class state, class action, class BloomFilter, bool verbose>
bool MBBDS<state, action, BloomFilter, verbose>::checkState(state midState)
{
	if (listReady) {
		if (middleStates.find(midState) != middleStates.end()) {
			return true;
		}
	}
	else if (firstRun || previousBloomfilter.contains(midState)) {
		if (outOfSpace) {
			currentBloomfilter.insert(midState);
		}
		else {
			if (middleStates.size() >= (int)(statesQuantityBound/2)) {
				outOfSpace = true;
				currentBloomfilter = BloomFilter(memoryBound/2, 1, previousBloomfilter.hashOffset + previousBloomfilter.getK());
				for (state possibleMidState : middleStates) {
					currentBloomfilter.insert(possibleMidState);
				}
				currentBloomfilter.insert(midState);
				middleStates.clear();
			}
			else {
				middleStates.insert(midState);
			}
		}
	}
	return false;
}

template<class state, class action, class BloomFilter, bool verbose>
void MBBDS<state, action, BloomFilter, verbose>::UpdateNextBound(double fCost)
{
	if (!fgreater(nextBound, fBound))
	{
		nextBound = fCost;
	}
	else if (fgreater(fCost, fBound) && fless(fCost, nextBound))
	{
		nextBound = fCost;
	}
}
#endif