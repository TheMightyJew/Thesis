/*
 *  IDMM.h
 *
 *  Created by Steven Danishevski on MARCH 20.
 *
 */

#ifndef IDMM_H
#define IDMM_H

#include <iostream>
#include <unordered_set>
#include "SearchEnvironment.h"

template <class state, class action, bool verbose = true>
class IDMM {
public:
	virtual ~IDMM() {}
	bool GetMidState(SearchEnvironment<state, action>* env, state fromState, state toState, state &midState, int secondsLimit=600);
	double getPathLength()	{ return pathLength; }
	uint64_t GetNodesExpanded() { return nodesExpanded; }
	uint64_t GetNodesTouched() { return nodesTouched; }
	void ResetNodeCount() { nodesExpanded = nodesTouched = 0; }
	unsigned long getDMMExpansions() { return dMMExpansions; }

private:
	unsigned long nodesExpanded, nodesTouched, dMMExpansions;
	double backwardBound;
	double forwardBound;
	double pathLength;
	bool DoIterationForward(SearchEnvironment<state, action>* env, state parent, state currState, double g, state& midState);	
	bool DoIterationBackward(SearchEnvironment<state, action>* env, state parent, state currState, double g, state& midState, state possibleMidState);
	
	state originGoal;
	state originStart;



};

template <class state, class action, bool verbose>
bool IDMM<state, action, verbose>::GetMidState(SearchEnvironment<state, action>* env,
	state fromState, state toState, state &midState, int secondsLimit)
{
	nodesExpanded = nodesTouched = 0;
	originStart = fromState;
	originGoal = toState;
	double heuristic = env->HCost(fromState, toState);
	backwardBound = (int)(heuristic / 2);
	forwardBound = heuristic - backwardBound;
	unsigned long nodesExpandedSoFar = 0;
	unsigned long lastIterationExpansions = 0;
	auto startTime = std::chrono::steady_clock::now();
	while (true){
		auto currentTime = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = currentTime-startTime;
		if(elapsed_seconds.count() >= secondsLimit){
			return false;
		}
		if (verbose){
			printf("\t\tBounds: %1.1f and %1.1f: ", forwardBound, backwardBound);
		}
		bool solved = DoIterationForward(env, originStart, originStart, 0, midState);
		lastIterationExpansions = nodesExpanded-nodesExpandedSoFar;
		if(verbose){
			printf("Nodes expanded: %d(%d)\n", lastIterationExpansions, nodesExpanded);
		}
		nodesExpandedSoFar = nodesExpanded;
		if (solved) {
			dMMExpansions = lastIterationExpansions;
			pathLength = backwardBound + forwardBound;
			return true;
		}
		else{
			if (forwardBound > backwardBound) {
				backwardBound = forwardBound;
			}
			else{
				forwardBound++;
			}			
		}
	}
	return false;
}


template <class state, class action, bool verbose>
bool IDMM<state, action, verbose>::DoIterationForward(SearchEnvironment<state, action>* env,
	state parent, state currState, double g, state& midState)
{
	double h = env->HCost(currState, originGoal);
	if (g > forwardBound  || fgreater(g + h, backwardBound + forwardBound))//not sure about that
	{
		return false;
	}
	else if (g == forwardBound) {
		if (DoIterationBackward(env, originGoal, originGoal, 0, midState, currState)) {
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
	for (unsigned int x = 0; x < neighbors.size(); x++)
	{
		if (neighbors[x] == parent) {
			continue;
		}
		double edgeCost = env->GCost(currState, neighbors[x]);
		if (DoIterationForward(env, currState, neighbors[x], g + edgeCost, midState)) {
			return true;
		}
	}
	return false;
}

template <class state, class action, bool verbose>
bool IDMM<state, action, verbose>::DoIterationBackward(SearchEnvironment<state, action>* env,
	state parent, state currState, double g, state& midState, state possibleMidState)
{
	double h = env->HCost(currState, possibleMidState);
	double originalH = env->HCost(currState, originStart);
	if (fgreater(g + h, backwardBound) || fgreater(g + originalH, forwardBound+backwardBound))//not sure about that
	{
		return false;
	}
	else if (currState == possibleMidState) {
		return true;
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
		if (DoIterationBackward(env, currState, neighbors[x], g + edgeCost, midState, possibleMidState)) {
			return true;
		}
	}
	return false;
}

#endif