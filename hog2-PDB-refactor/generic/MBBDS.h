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
#include "PancakeHasher.h"
#include <math.h>

template <class state, class action, class BloomFilter, bool verbose = true>
class MBBDS {
public:
	MBBDS(unsigned long statesQuantityBound) {
		this->statesQuantityBound = statesQuantityBound;
	}
	virtual ~MBBDS() {}
	bool GetMidState(SearchEnvironment<state, action>* env, state fromState, state toState, state &midState, int secondsLimit=600, double startingFBound=0);
	double getPathLength()	{ return pathLength; }
	uint64_t GetNodesExpanded() { return nodesExpanded; }
	uint64_t GetNodesTouched() { return nodesTouched; }
	double getLastBound() { return backwardBound+forwardBound; }
	void ResetNodeCount() { nodesExpanded = nodesTouched = 0; }
private:
	unsigned long nodesExpanded, nodesTouched, statesQuantityBound;
	BloomFilter previousBloomfilter;
	BloomFilter currentBloomfilter;
	bool listReady;
	bool outOfSpace;
	bool firstRun;
	double backwardBound;
	double forwardBound;
	double pathLength;
	std::unordered_set<state> middleStates;
	unsigned long memoryBound;

	bool DoIteration(SearchEnvironment<state, action>* env,
		state parent, state currState, double bound, double g, state& midState);
	bool checkState(state midState);
	
	state goal;
	state start;



};

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
	state from;
	double initialHeuristic = env->HCost(fromState, toState);
	startingFBound = round(std::max(initialHeuristic, startingFBound));
	backwardBound = (int)(startingFBound / 2);
	forwardBound = startingFBound - backwardBound;
	bool forwardSearch;
	int saturationIncreased = 0;
	int saturationMaxIncreasements = 100;
	int iteration_num = 0;
	double bound;
	firstRun = true;
	forwardSearch = true;
	listReady = false;
	middleStates.clear();
	unsigned long nodesExpandedSoFar = 0;
	double last_saturation = 1;
	auto startTime = std::chrono::steady_clock::now();
	while(true){
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
				bound = forwardBound;
				from = fromState;
				goal = toState;
			}
			else {
				bound = backwardBound;
				from = toState;
				goal = fromState;
			}
			outOfSpace = false;
			bool solved = DoIteration(env, from, from, bound, 0, midState);
			if(verbose){
				printf("Nodes expanded: %d(%d)\n", nodesExpanded-nodesExpandedSoFar, nodesExpanded);
			}
			nodesExpandedSoFar = nodesExpanded;
			iteration_num++;
			if (solved) {
				pathLength = backwardBound + forwardBound;
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
		}	
		
		last_saturation = 1;
		saturationIncreased = 0;
		previousBloomfilter.clear();
		firstRun = true;
		forwardSearch = true;
		listReady = false;
		middleStates.clear();
		if (forwardBound == backwardBound) {
			forwardBound++;
		}
		else{
			forwardBound = backwardBound = std::max(forwardBound, backwardBound);
		}
	}
	return false;
}


template <class state, class action, class BloomFilter, bool verbose>
bool MBBDS<state, action, BloomFilter, verbose>::DoIteration(SearchEnvironment<state, action>* env,
	state parent, state currState, double bound, double g, state& midState)
{
	double h = env->HCost(currState, goal);
	if (g == bound) {
		bool solution = checkState(currState);
		if (solution) {
			midState = currState;
			return true;
		}
		else{
			return false;
		}
	}
	else if (g > bound  || fgreater(g + h, backwardBound + forwardBound))//not sure about that
	{
		return false;
	}
	std::vector<state> neighbors;
	env->GetSuccessors(currState, neighbors);
	nodesTouched += neighbors.size();
	nodesExpanded++;
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
			if (middleStates.size() == (int)(statesQuantityBound/2)) {
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

#endif