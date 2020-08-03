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
#include "AStarOpenClosed.h"
#include "MM.h"
#include <math.h>



template <class state, class action, bool verbose = true>
class IDMM {
public:
	IDMM(bool front2frontH=false) { this->front2frontH = front2frontH;}
	virtual ~IDMM() {}
	bool GetMidState(SearchEnvironment<state, action>* env, state fromState, state toState, state &midState, int secondsLimit=600, double startingFBound=0);
	bool GetMidStateFromLists(SearchEnvironment<state, action>* env, state fromState, state toState, state &midState, int secondsLimit=600, double startingFBound=0, AStarOpenClosed<state, MMCompare<state>> forwardList = AStarOpenClosed<state, MMCompare<state>>(), AStarOpenClosed<state, MMCompare<state>> backwardList = AStarOpenClosed<state, MMCompare<state>>());
	bool GetMidStateFromForwardList(SearchEnvironment<state, action>* env, state fromState, state toState, state &midState, int secondsLimit=600, AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> forwardList = AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>>());
	double getPathLength()	{ return pathLength; }
	uint64_t GetNodesExpanded() { return nodesExpanded; }
	uint64_t GetNecessaryExpansions() { return necessaryExpansions; }
	uint64_t GetNodesTouched() { return nodesTouched; }
	void ResetNodeCount() { nodesExpanded = nodesTouched = 0; }
	unsigned long getDMMExpansions() { return dMMExpansions; }

private:
	unsigned long nodesExpanded, nodesTouched, dMMExpansions;
	double backwardBound;
	double forwardBound;
	double pathLength = std::numeric_limits<double>::max();
	bool DoIterationForward(SearchEnvironment<state, action>* env, state parent, state currState, double g, state& midState);	
	bool DoIterationBackward(SearchEnvironment<state, action>* env, state parent, state currState, double g, state& midState, state possibleMidState);
	void updateBoundsG(double minOpenG, double maxOpenG);
	void buildMatrix(std::vector<AStarOpenClosedData<state>> &openList, std::vector<std::vector<AStarOpenClosedData<state>>> &matrix);
	
	state originGoal;
	state originStart;
	unsigned long dMMLastIterExpansions = 0;
	unsigned long necessaryExpansions = 0;
	
	AStarOpenClosed<state, MMCompare<state>, AStarOpenClosedData<state>> forwardList;
	AStarOpenClosed<state, MMCompare<state>, AStarOpenClosedData<state>> backwardList;
	std::vector<AStarOpenClosedData<state>> forwardOpenList;
	std::vector<AStarOpenClosedData<state>> backwardOpenList;
	std::vector<std::vector<AStarOpenClosedData<state>>> forwardMatrix;
	std::vector<std::vector<AStarOpenClosedData<state>>> backwardMatrix;
	bool readyOpenLists;
	bool front2frontH;
	bool firstBounds;


};
template <class state, class action, bool verbose>
bool IDMM<state, action, verbose>::GetMidStateFromForwardList(SearchEnvironment<state, action>* env,
	state fromState, state toState, state &midState, int secondsLimit, AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> forwardList)
{
	AStarOpenClosed<state, MMCompare<state>, AStarOpenClosedData<state>> newForwardList;
	double minF = std::numeric_limits<double>::max();
	for (AStarOpenClosedDataWithF<state> forwardState : forwardList.getElements()){
		if(forwardState.where == kOpenList){
			minF = std::min(minF, forwardState.h+forwardState.g);
			newForwardList.AddOpenNode(forwardState.data, env->GetStateHash(forwardState.data), forwardState.g, forwardState.h);
		}
		else if(forwardState.where == kClosedList){
			newForwardList.AddClosedNode(forwardState.data, env->GetStateHash(forwardState.data), forwardState.g, forwardState.h);
		}
	}
	
	AStarOpenClosed<state, MMCompare<state>, AStarOpenClosedData<state>> newBackwardList;
	newBackwardList.AddOpenNode(toState, env->GetStateHash(toState), 0, env->HCost(toState, fromState));
	
	return GetMidStateFromLists(env, fromState, toState, midState, secondsLimit, minF, newForwardList, newBackwardList);
}

template <class state, class action, bool verbose>
void IDMM<state, action, verbose>::buildMatrix(std::vector<AStarOpenClosedData<state>> &openList, std::vector<std::vector<AStarOpenClosedData<state>>> &matrix){
	sort( openList.begin( ), openList.end( ), [ ]( const AStarOpenClosedData<state>& lhs, const AStarOpenClosedData<state>& rhs )
	{
	   return lhs.h+lhs.g < rhs.h+rhs.g;
	});
	
	double f = -1;
	for(AStarOpenClosedData<state> openState:openList){
		if(openState.h+openState.g > f){
			f = openState.h+openState.g;
			std::vector<AStarOpenClosedData<state>> newVector;
			matrix.push_back(newVector);
		}
		matrix.back().push_back(openState);
	}
	for(std::vector<int>::size_type i = 0; i != matrix.size(); i++) {
		sort( matrix[i].begin( ), matrix[i].end( ), [ ]( const AStarOpenClosedData<state>& lhs, const AStarOpenClosedData<state>& rhs )
		{
		   return lhs.g < rhs.g;
		});
	}
}

template <class state, class action, bool verbose>
bool IDMM<state, action, verbose>::GetMidStateFromLists(SearchEnvironment<state, action>* env,
	state fromState, state toState, state &midState, int secondsLimit, double startingFBound, AStarOpenClosed<state, MMCompare<state>, AStarOpenClosedData<state>> forwardList, AStarOpenClosed<state, MMCompare<state>, AStarOpenClosedData<state>> backwardList)
{
	auto startTime = std::chrono::steady_clock::now();
	this->readyOpenLists = true;
	this->forwardList = forwardList;
	this->backwardList = backwardList;
	nodesExpanded = nodesTouched = 0;
	originStart = fromState;
	originGoal = toState;
	double minF = std::numeric_limits<double>::max();
	double maxOpenG = 0;
	double minOpenG = std::numeric_limits<double>::max();

	for (int x = 0; x < forwardList.OpenSize(); x++){
		AStarOpenClosedData<state> openState = forwardList.getElements()[forwardList.GetOpenItem(x)];
		
		minOpenG = std::min(minOpenG, openState.g);
		maxOpenG = std::max(maxOpenG , openState.g);
			
		forwardOpenList.push_back(openState);
	}
	buildMatrix(forwardOpenList, forwardMatrix);
	minF = forwardOpenList.front().g + forwardOpenList.front().h;

	
	for (int x = 0; x < backwardList.OpenSize(); x++){
		AStarOpenClosedData<state> openState = backwardList.getElements()[backwardList.GetOpenItem(x)];
		backwardOpenList.push_back(openState);
	}
	buildMatrix(backwardOpenList, backwardMatrix);

	double initialHeuristic = env->HCost(fromState, toState);
	startingFBound = std::max(startingFBound, std::max(initialHeuristic, minF));
	forwardBound = std::max(ceil(startingFBound / 2), minOpenG);
	backwardBound = startingFBound - forwardBound;
	firstBounds = true;
	unsigned long nodesExpandedSoFar = 0;
	unsigned long previousIterationExpansions = 0;
	while (true){
		dMMLastIterExpansions = 0;
		auto currentTime = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = currentTime-startTime;
		if(elapsed_seconds.count() >= secondsLimit){
			return false;
		}
		if (verbose){
			printf("\t\tBounds: %1.1f and %1.1f: ", forwardBound, backwardBound);
		}
		bool solved = false;
		for(std::vector<AStarOpenClosedData<state>> fVector:forwardMatrix){
			if(fVector.front().g+fVector.front().h > forwardBound+backwardBound){
				break;
			}
			for (AStarOpenClosedData<state> openState: fVector){
				auto currentTime = std::chrono::steady_clock::now();
				std::chrono::duration<double> elapsed_seconds = currentTime-startTime;
				if(elapsed_seconds.count() >= secondsLimit){
					return false;
				}
				if(openState.g > forwardBound)
					break;
				if(firstBounds || openState.g == forwardBound){
					solved = DoIterationForward(env, openState.data, openState.data, openState.g, midState);
				}
				if(solved){
					break;
				}
			}
			if(solved){
				break;
			}
		}
		if(verbose){
			printf("Nodes expanded: %d(%d)\n", nodesExpanded-nodesExpandedSoFar, nodesExpanded);
		}
		if (solved) {
			dMMExpansions = previousIterationExpansions + dMMLastIterExpansions;
			necessaryExpansions += nodesExpandedSoFar;
			pathLength = std::min(pathLength, backwardBound + forwardBound);
			return true;
		}
		else{
			updateBoundsG(minOpenG, maxOpenG);			
		}
		previousIterationExpansions = nodesExpanded-nodesExpandedSoFar;
		nodesExpandedSoFar = nodesExpanded;
	}
	return false;
}

template <class state, class action, bool verbose>
void IDMM<state, action, verbose>::updateBoundsG(double minOpenG, double maxOpenG){
	/*if(forwardBound>backwardBound){
		backwardBound++;
	}
	else{
		forwardBound++;
	}*/
	if(maxOpenG > forwardBound && backwardBound-1 >= 0){
		forwardBound++;
		backwardBound--;
		firstBounds = false;
	}
	else{
		double fullBound = forwardBound + backwardBound + 1;
		forwardBound = std::max(minOpenG, ceil(fullBound/2));
		backwardBound = fullBound - forwardBound;
		firstBounds = true;
	}
}
template <class state, class action, bool verbose>
bool IDMM<state, action, verbose>::GetMidState(SearchEnvironment<state, action>* env,
	state fromState, state toState, state &midState, int secondsLimit, double startingFBound)
{
	auto startTime = std::chrono::steady_clock::now();
	this->readyOpenLists = false;
	this->forwardList = forwardList;
	this->backwardList = backwardList;
	nodesExpanded = nodesTouched = 0;
	originStart = fromState;
	originGoal = toState;
	double initialHeuristic = env->HCost(fromState, toState);
	startingFBound = std::max(initialHeuristic, startingFBound); 
	forwardBound = ceil(startingFBound / 2);
	backwardBound = startingFBound - forwardBound;
	unsigned long nodesExpandedSoFar = 0;
	unsigned long previousIterationExpansions = 0;
	while (true){
		dMMLastIterExpansions = 0;
		auto currentTime = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = currentTime-startTime;
		if(elapsed_seconds.count() >= secondsLimit){
			return false;
		}
		if (verbose){
			printf("\t\tBounds: %1.1f and %1.1f: ", forwardBound, backwardBound);
		}
		bool solved = DoIterationForward(env, originStart, originStart, 0, midState);
		if(verbose){
			printf("Nodes expanded: %d(%d)\n", nodesExpanded-nodesExpandedSoFar, nodesExpanded);
		}
		if (solved) {
			dMMExpansions = previousIterationExpansions + dMMLastIterExpansions;
			necessaryExpansions += nodesExpandedSoFar;
			pathLength = std::min(pathLength, backwardBound + forwardBound);
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
		previousIterationExpansions = nodesExpanded-nodesExpandedSoFar;
		nodesExpandedSoFar = nodesExpanded;
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
		if(readyOpenLists){	
			for(std::vector<AStarOpenClosedData<state>> fVector:backwardMatrix){
				if(fVector.front().g+fVector.front().h > forwardBound+backwardBound){
					break;
				}
				for (AStarOpenClosedData<state> openState: fVector){
					if(openState.g > backwardBound)
						break;
					if (DoIterationBackward(env, openState.data, openState.data, openState.g, midState, currState)) {
						pathLength += g;
						midState = currState;
						return true;
					}
				}
			}
			return false;
		}
		else if (DoIterationBackward(env, originGoal, originGoal, 0, midState, currState)) {
			pathLength += g;
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
	if(g + h == forwardBound+backwardBound){
		dMMLastIterExpansions++;
	}
	for (unsigned int x = 0; x < neighbors.size(); x++)
	{
		uint64_t childID;
		if (neighbors[x] == parent || forwardList.Lookup(env->GetStateHash(neighbors[x]), childID) == kClosedList) {
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
	if ((front2frontH && fgreater(g + h, backwardBound)) || g>backwardBound || fgreater(g + originalH, forwardBound+backwardBound)){
		return false;
	}
	else if (currState == possibleMidState) {
		pathLength = g;
		return true;
	}
	std::vector<state> neighbors;
	env->GetSuccessors(currState, neighbors);
	nodesTouched += neighbors.size();
	nodesExpanded++;
	if(g + originalH == forwardBound+backwardBound){
		dMMLastIterExpansions++;
	}
	for (unsigned int x = 0; x < neighbors.size(); x++)
	{
		uint64_t childID;
		if (neighbors[x] == parent || backwardList.Lookup(env->GetStateHash(neighbors[x]), childID) == kClosedList) {
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