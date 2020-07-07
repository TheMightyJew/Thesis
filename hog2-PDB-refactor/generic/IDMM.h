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



template <class state, class action, bool verbose = true>
class IDMM {
public:
	IDMM(bool front2frontH=false) { this->front2frontH = front2frontH;}
	virtual ~IDMM() {}
	bool GetMidState(SearchEnvironment<state, action>* env, state fromState, state toState, state &midState, int secondsLimit=600, double startingFBound=0, bool readyOpenLists=false, AStarOpenClosed<state, MMCompare<state>> forwardList = AStarOpenClosed<state, MMCompare<state>>(), AStarOpenClosed<state, MMCompare<state>> backwardList = AStarOpenClosed<state, MMCompare<state>>());
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
	double pathLength;
	bool DoIterationForward(SearchEnvironment<state, action>* env, state parent, state currState, double g, state& midState);	
	bool DoIterationBackward(SearchEnvironment<state, action>* env, state parent, state currState, double g, state& midState, state possibleMidState);
	
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


};

template <class state, class action, bool verbose>
bool IDMM<state, action, verbose>::GetMidState(SearchEnvironment<state, action>* env,
	state fromState, state toState, state &midState, int secondsLimit, double startingFBound, bool readyOpenLists, AStarOpenClosed<state, MMCompare<state>, AStarOpenClosedData<state>> forwardList, AStarOpenClosed<state, MMCompare<state>, AStarOpenClosedData<state>> backwardList)
{
	auto startTime = std::chrono::steady_clock::now();
	this->readyOpenLists = readyOpenLists;
	this->forwardList = forwardList;
	this->backwardList = backwardList;
	nodesExpanded = nodesTouched = 0;
	originStart = fromState;
	originGoal = toState;
	if(readyOpenLists){
		double minF = 0;

		for (int x = 0; x < forwardList.OpenSize(); x++){
			AStarOpenClosedData<state> openState = forwardList.getElements()[forwardList.GetOpenItem(x)];
			double fValue = openState.h + openState.g;
			if(minF == 0 || minF > fValue){
					minF = fValue;
			}
		}
		startingFBound = std::max(minF, startingFBound);

		for (int x = 0; x < forwardList.OpenSize(); x++){
			AStarOpenClosedData<state> openState = forwardList.getElements()[forwardList.GetOpenItem(x)];
			forwardOpenList.push_back(openState);
		}
		sort( forwardOpenList.begin( ), forwardOpenList.end( ), [ ]( const AStarOpenClosedData<state>& lhs, const AStarOpenClosedData<state>& rhs )
		{
		   return lhs.h+lhs.g < rhs.h+rhs.g;
		});

		for (int x = 0; x < backwardList.OpenSize(); x++){
			AStarOpenClosedData<state> openState = backwardList.getElements()[backwardList.GetOpenItem(x)];
			backwardOpenList.push_back(openState);
		}
		sort( backwardOpenList.begin( ), backwardOpenList.end( ), [ ]( const AStarOpenClosedData<state>& lhs, const AStarOpenClosedData<state>& rhs )
		{
		   return lhs.h+lhs.g < rhs.h+rhs.g;
		});
		double f = -1;
		for(AStarOpenClosedData<state> openState:forwardOpenList){
			if(openState.h+openState.g > f){
				f = openState.h+openState.g;
				std::vector<AStarOpenClosedData<state>> newVector;
				forwardMatrix.push_back(newVector);
			}
			forwardMatrix.back().push_back(openState);
		}
		for(std::vector<int>::size_type i = 0; i != forwardMatrix.size(); i++) {
			sort( forwardMatrix[i].begin( ), forwardMatrix[i].end( ), [ ]( const AStarOpenClosedData<state>& lhs, const AStarOpenClosedData<state>& rhs )
			{
			   return lhs.g < rhs.g;
			});
		}
	
		f = -1;
		for(AStarOpenClosedData<state> openState:backwardOpenList){
			if(openState.h+openState.g > f){
				f = openState.h+openState.g;
				std::vector<AStarOpenClosedData<state>> newVector;
				backwardMatrix.push_back(newVector);
			}
			backwardMatrix.back().push_back(openState);
		}
		for(std::vector<int>::size_type i = 0; i != backwardMatrix.size(); i++) {
			sort( backwardMatrix[i].begin( ), backwardMatrix[i].end( ), [ ]( const AStarOpenClosedData<state>& lhs, const AStarOpenClosedData<state>& rhs )
			{
			   return lhs.g < rhs.g;
			});
		}

		/*auto currentTime = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = currentTime-startTime;
		std::cout << elapsed_seconds.count() << std::endl;*/
		//necessaryExpansions = 
	}
	double initialHeuristic = env->HCost(fromState, toState);
	startingFBound = std::max(initialHeuristic, startingFBound); 
	backwardBound = (int)(startingFBound / 2);
	forwardBound = startingFBound - backwardBound;
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
		if(readyOpenLists){
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
					solved = DoIterationForward(env, openState.data, openState.data, openState.g, midState);
					if(solved){
						break;
					}
				}
				if(solved){
					break;
				}
			}
			/*for (AStarOpenClosedData<state> openState: forwardOpenList){
				if(openState.g+openState.h > forwardBound+backwardBound)
					break;
				solved = DoIterationForward(env, openState.data, openState.data, openState.g, midState);
				if(solved){
					break;
				}
			}*/
		}
		else{
			solved = DoIterationForward(env, originStart, originStart, 0, midState);
		}
		if(verbose){
			printf("Nodes expanded: %d(%d)\n", nodesExpanded-nodesExpandedSoFar, nodesExpanded);
		}
		if (solved) {
			dMMExpansions = previousIterationExpansions + dMMLastIterExpansions;
			necessaryExpansions += nodesExpandedSoFar;
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
			/*for (AStarOpenClosedData<state> openState: backwardOpenList){
				if(openState.g+openState.h > forwardBound+backwardBound)
					break;
				if (DoIterationBackward(env, openState.data, openState.data, openState.g, midState, currState)) {
					midState = currState;
					return true;
				}
			}
			return false;*/
			
			for(std::vector<AStarOpenClosedData<state>> fVector:backwardMatrix){
				if(fVector.front().g+fVector.front().h > forwardBound+backwardBound){
					break;
				}
				for (AStarOpenClosedData<state> openState: fVector){
					if(openState.g > backwardBound)
						break;
					if (DoIterationBackward(env, openState.data, openState.data, openState.g, midState, currState)) {
						midState = currState;
						return true;
					}
				}
			}
			return false;
		}
		else if (DoIterationBackward(env, originGoal, originGoal, 0, midState, currState)) {
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
	if ((front2frontH && fgreater(g + h, backwardBound)) || g>backwardBound || fgreater(g + originalH, forwardBound+backwardBound))//not sure about that
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