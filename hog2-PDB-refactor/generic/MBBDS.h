/*
 *  MBBDS.h
 *
 *  Created by Steven Danishevski on DEC 19.
 *
 */

#ifndef MBBDS_H
#define MBBDS_H

#include <iostream>
#include <unordered_set>
#include "SearchEnvironment.h"
#include "PancakeHasher.h"


template <class state, class action, class BloomFilter, bool verbose = true>
class MBBDS {
public:
	MBBDS(unsigned long statesQuantityBound, long long unsigned int memoryBound) {
		this->statesQuantityBound = statesQuantityBound;
		this->memoryBound = memoryBound;
	}
	virtual ~MBBDS() {}
	double GetMidState(SearchEnvironment<state, action>* env, state fromState, state toState, state &midState);
	bool checkState(state midState);
	uint64_t GetNodesExpanded() { return nodesExpanded; }
	uint64_t GetNodesTouched() { return nodesTouched; }
	void ResetNodeCount() { nodesExpanded = nodesTouched = 0; }
private:
	unsigned long nodesExpanded, nodesTouched, statesQuantityBound;
	long long unsigned int memoryBound;
	BloomFilter previousBloomfilter;
	BloomFilter currentBloomfilter;
	bool listReady;
	bool outOfSpace;
	bool firstRun;
	double backwardBound;
	double forwardBound;
	std::unordered_set<state> middleStates;

	bool DoIteration(SearchEnvironment<state, action>* env,
		state parent, state currState, double bound, double g, state& midState);

	state goal;
	state start;



};

template <class state, class action, class BloomFilter, bool verbose>
double MBBDS<state, action, BloomFilter, verbose>::GetMidState(SearchEnvironment<state, action>* env,
	state fromState, state toState, state &midState)
{
	if(verbose){
		printf("\nStarting to solve with MBBDS\n");
	}
	nodesExpanded = nodesTouched = 0;
	state from;
	double heuristic = env->HCost(fromState, toState);
	backwardBound = (int)(heuristic / 2);
	forwardBound = heuristic - backwardBound;
	bool forwardSearch;
	int iteration_num = 0;
	double bound;
	firstRun = true;
	forwardSearch = true;
	listReady = false;
	middleStates.clear();
	unsigned long nodesExpandedSoFar = 0;
	while(true){
		if (verbose){
			printf("Bounds: %f and %f\n", forwardBound, backwardBound);
		}
		while (true){
			if (verbose){
				printf("Starting iteration number: %d. ", iteration_num);
				if(false){
					printf("start middleStates.size() = %d\n", middleStates.size());
					printf("start previousBloomfilter.getSaturation() = %f\n", previousBloomfilter.getSaturation());
					printf("start listReady = %d\n", listReady);
					printf("start firstRun = %d\n", firstRun);
					printf("start forwardSearch = %d\n", forwardSearch);
				}
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
			//env->StoreGoal(goal);
			outOfSpace = false;
			bool solved = DoIteration(env, from, from, bound, 0, midState);
			if(verbose){
				printf("Nodes expanded: %d(%d)\n", nodesExpanded-nodesExpandedSoFar, nodesExpanded);
			}
			nodesExpandedSoFar = nodesExpanded;
			iteration_num++;
			if (solved) {
				return backwardBound + forwardBound;
			}
			if (listReady) {// no solution found
				if(verbose){
					std::cout << "No solution" << std::endl;
				}
				break;
			}
			if (!outOfSpace) {
				listReady = true;
				previousBloomfilter.clear();
			}
			else {
				double saturation = currentBloomfilter.getSaturation();
				if (saturation == 1) {
					if(verbose){
						std::cout << "fuck" << std::endl;
					}
					return -1; //bloomfilter is fluded.
				}
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
		
		if (forwardBound > backwardBound) {
			backwardBound = forwardBound;
		}
		else{
			forwardBound++;
		}
		firstRun = true;
		forwardSearch = true;
		listReady = false;
		middleStates.clear();
		previousBloomfilter.clear();
	}
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
		if (DoIteration(env, currState, neighbors[x], bound, g + edgeCost,midState)) {
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
			if (middleStates.size() == statesQuantityBound) {
				outOfSpace = true;
				currentBloomfilter = BloomFilter(previousBloomfilter.hashOffset + previousBloomfilter.getK());
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