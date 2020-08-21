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
	IDMM(bool front2frontH=false, bool isConsistent = false, double smallestEdge=1) { this->front2frontH = front2frontH; this->isConsistent = isConsistent; this->smallestEdge = smallestEdge;}
	virtual ~IDMM() {}
	bool GetMidState(SearchEnvironment<state, action>* env, state fromState, state toState, state &midState, int secondsLimit=600, double startingFBound=0);
	bool GetMidStateFromLists(SearchEnvironment<state, action>* env, state fromState, state toState, state &midState, int secondsLimit=600, double startingFBound=0, AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> forwardList = AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>>(), AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> backwardList = AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>>(), bool detectDuplicates = true);
	bool GetMidStateFromForwardList(SearchEnvironment<state, action>* env, state fromState, state toState, state &midState, int secondsLimit=600, AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> forwardList = AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>>(),bool detectDuplicates = true);
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
	double fBound;
	double pathLength = std::numeric_limits<double>::max();
  bool detectDuplicates;
	bool DoIterationForward(SearchEnvironment<state, action>* env, state parent, state currState, double g, state& midState);	
	bool DoIterationBackward(SearchEnvironment<state, action>* env, state parent, state currState, double g, state& midState, state possibleMidState, double possibleMidStateG, double otherH = 0,double otherError = 0);
	void buildMatrix(std::vector<AStarOpenClosedDataWithF<state>> &openList, std::vector<std::vector<AStarOpenClosedDataWithF<state>>> &matrix);
	double nextBound[2];
  int const forwardLoc = 0;
  int const backwardLoc = 0;
	void UpdateNextBound(double currBound, double fCost, int loc);
	state originGoal;
	state originStart;
	unsigned long dMMLastIterExpansions = 0;
	unsigned long necessaryExpansions = 0;
	unsigned long nodesExpandedSoFar = 0;
	unsigned long previousIterationExpansions = 0;
	
	AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> forwardList;
	AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> backwardList;
	std::vector<AStarOpenClosedDataWithF<state>> forwardOpenList;
	std::vector<AStarOpenClosedDataWithF<state>> backwardOpenList;
	std::vector<std::vector<AStarOpenClosedDataWithF<state>>> forwardMatrix;
	std::vector<std::vector<AStarOpenClosedDataWithF<state>>> backwardMatrix;
	bool readyOpenLists;
	bool front2frontH;
	bool firstBounds;
	bool isConsistent;
	double minForwardError = 0;
	double minBackwardError = 0;
	double smallestEdge;
	std::chrono::steady_clock::time_point startTimeTest;
	//double prevMinForwardError = 0;
	//double prevMinBackwardError = 0;

};
template <class state, class action, bool verbose>
bool IDMM<state, action, verbose>::GetMidStateFromForwardList(SearchEnvironment<state, action>* env,
	state fromState, state toState, state &midState, int secondsLimit, AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> forwardList, bool detectDuplicates)
{
	startTimeTest = std::chrono::steady_clock::now();
	
	AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> newBackwardList;
	double g = 0;
	double h = env->HCost(toState, fromState);
	newBackwardList.AddOpenNode(toState, env->GetStateHash(toState), g+h, 0, h);
	
	return GetMidStateFromLists(env, fromState, toState, midState, secondsLimit, 0, forwardList, newBackwardList, detectDuplicates);
}

template <class state, class action, bool verbose>
void IDMM<state, action, verbose>::buildMatrix(std::vector<AStarOpenClosedDataWithF<state>> &openList, std::vector<std::vector<AStarOpenClosedDataWithF<state>>> &matrix){
	sort( openList.begin( ), openList.end( ), [ ]( const AStarOpenClosedDataWithF<state>& lhs, const AStarOpenClosedDataWithF<state>& rhs )
	{
	   return lhs.h+lhs.g < rhs.h+rhs.g;
	});
	
	double f = -1;
	for(AStarOpenClosedDataWithF<state> openState:openList){
		if(openState.h+openState.g > f){
			f = openState.h+openState.g;
			std::vector<AStarOpenClosedDataWithF<state>> newVector;
			matrix.push_back(newVector);
		}
		matrix.back().push_back(openState);
	}
	for(std::vector<int>::size_type i = 0; i != matrix.size(); i++) {
		sort( matrix[i].begin( ), matrix[i].end( ), [ ]( const AStarOpenClosedDataWithF<state>& lhs, const AStarOpenClosedDataWithF<state>& rhs )
		{
		   return lhs.g < rhs.g;
		});
	}
}

template <class state, class action, bool verbose>
bool IDMM<state, action, verbose>::GetMidStateFromLists(SearchEnvironment<state, action>* env,
	state fromState, state toState, state &midState, int secondsLimit, double startingFBound, AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> forwardList, AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> backwardList, bool detectDuplicates)
{
	auto startTime = std::chrono::steady_clock::now();
	this->readyOpenLists = true;
	this->forwardList = forwardList;
	this->backwardList = backwardList;
	this->detectDuplicates = detectDuplicates;
	nodesExpanded = nodesTouched = 0;
	originStart = fromState;
	originGoal = toState;
	double minF = std::numeric_limits<double>::max();
	double minOpenG = std::numeric_limits<double>::max();
	minForwardError = DBL_MAX;
	minBackwardError = DBL_MAX;
  
	for (int x = 0; x < forwardList.OpenSize(); x++){
		AStarOpenClosedDataWithF<state> openState = forwardList.getElements()[forwardList.GetOpenItem(x)];
		minOpenG = std::min(minOpenG, openState.g);
		minF = std::min(minF, openState.f);
		minForwardError = std::min(minForwardError, openState.g - env->HCost(originStart, openState.data));
		forwardOpenList.push_back(openState);
	}

	for (int x = 0; x < backwardList.OpenSize(); x++){
		AStarOpenClosedDataWithF<state> openState = backwardList.getElements()[backwardList.GetOpenItem(x)];
		minF = std::min(minF, openState.f);
		minBackwardError = std::min(minBackwardError, openState.g - env->HCost(openState.data,originGoal));
		backwardOpenList.push_back(openState);
	}
	
	if (!isConsistent){
		minBackwardError = minForwardError = 0;
	}

  	sort( forwardOpenList.begin( ), forwardOpenList.end( ), [ ]( const AStarOpenClosedDataWithF<state>& lhs, const AStarOpenClosedDataWithF<state>& rhs )
	{
	   return lhs.h < rhs.h;
	});
	sort( backwardOpenList.begin( ), backwardOpenList.end( ), [ ]( const AStarOpenClosedDataWithF<state>& lhs, const AStarOpenClosedDataWithF<state>& rhs )
	{
	   return lhs.h < rhs.h;
	});
	
	double initialHeuristic = env->HCost(fromState, toState);
	fBound = nextBound[forwardLoc] = nextBound[backwardLoc] = std::max(smallestEdge,std::max(startingFBound, std::max(initialHeuristic, minF)));
	forwardBound = std::max(minOpenG, ceil(fBound/2) - smallestEdge);
	firstBounds = true;
	auto curtime = std::chrono::steady_clock::now();
	std::chrono::duration<double> process_time = curtime-startTimeTest;
	printf("%f\n", process_time);
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
		for (AStarOpenClosedDataWithF<state> openState: forwardOpenList){
			auto currentTime = std::chrono::steady_clock::now();
			std::chrono::duration<double> elapsed_seconds = currentTime-startTime;
			if(elapsed_seconds.count() >= secondsLimit){
				return false;
			}
			solved = DoIterationForward(env, openState.data, openState.data, openState.g, midState);
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
			for (AStarOpenClosedDataWithF<state> forwardState : forwardList.getElements()){
				if(forwardState.where == kClosedList && forwardState.f < pathLength){
					necessaryExpansions++;
				}
			}
			for (AStarOpenClosedDataWithF<state> backwardState : backwardList.getElements()){
				if(backwardState.where == kClosedList && backwardState.f < pathLength){
					necessaryExpansions++;
				}
			}
			return true;
		}
		else{
			fBound = std::max(nextBound[forwardLoc],nextBound[backwardLoc]);
			forwardBound = std::max(minOpenG, ceil(fBound/2) - smallestEdge);
		}
	}
	return false;
}

template <class state, class action, bool verbose>
bool IDMM<state, action, verbose>::GetMidState(SearchEnvironment<state, action>* env,
	state fromState, state toState, state &midState, int secondsLimit, double startingFBound)
{
	auto startTime = std::chrono::steady_clock::now();
	this->readyOpenLists = false;
	nodesExpanded = nodesTouched = 0;
	originStart = fromState;
	originGoal = toState;
	double initialHeuristic = env->HCost(fromState, toState);
	fBound = 	nextBound[forwardLoc] = nextBound[backwardLoc] = std::max(smallestEdge,std::max(startingFBound, initialHeuristic));
	forwardBound = ceil(fBound / 2) - smallestEdge;
	//printf("%f|%f|%f\n", fBound, smallestEdge, forwardBound);
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
			return true;
		}
		else{
			fBound = std::max(nextBound[forwardLoc],nextBound[backwardLoc]);
			forwardBound = ceil(fBound / 2) - smallestEdge;
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
	if (fgreater(g + h + minBackwardError, fBound)){
		UpdateNextBound(fBound, g + h + minBackwardError,forwardLoc);
		return false;
	}
	else if (g > forwardBound) {
		backwardBound = fBound - g - smallestEdge;
		double error = 0;
		if (isConsistent){
		  error = g - env->HCost(currState,originStart);
		}
		if(readyOpenLists){	
			for (AStarOpenClosedDataWithF<state> openState: backwardOpenList){
				//if(openState.g + openState.h > fBound || openState.g > backwardBound){
				if(openState.g + openState.h > fBound){
				  UpdateNextBound(fBound, openState.g + openState.h,forwardLoc);
				  continue;
				}
				if (DoIterationBackward(env, openState.data, openState.data, openState.g, midState, currState, g, h, error)) {
				  pathLength += g;
				  midState = currState;
				  return true;
				}
			}
			return false;
		}
		else{
			if (DoIterationBackward(env, originGoal, originGoal, 0, midState, currState, g, h, error)) {
				pathLength += g;
				midState = currState;
				return true;
			}
			else{
				return false;
			}
		}
	}
	std::vector<state> neighbors;
	env->GetSuccessors(currState, neighbors);
	nodesTouched += neighbors.size();
	nodesExpanded++;
	if(g + h == fBound){
		dMMLastIterExpansions++;
	}
	for (unsigned int x = 0; x < neighbors.size(); x++){
		uint64_t childID;
		double edgeCost = env->GCost(currState, neighbors[x]);
		if (neighbors[x] == parent || (detectDuplicates && forwardList.getElements().size()>1 && forwardList.Lookup(env->GetStateHash(neighbors[x]), childID) != kNotFound && (forwardList.Lookup(childID).where == kClosedList || (forwardList.Lookup(childID).where == kOpenList && forwardList.Lookup(childID).g <= g + edgeCost)))) {
			continue;
		}
		if (DoIterationForward(env, currState, neighbors[x], g + edgeCost, midState)) {
			return true;
		}
	}
	return false;
}

template <class state, class action, bool verbose>
bool IDMM<state, action, verbose>::DoIterationBackward(SearchEnvironment<state, action>* env,
	state parent, state currState, double g, state& midState, state possibleMidState, double possibleMidStateG, double otherH, double otherError)
{	
	if (currState == possibleMidState && g + possibleMidStateG <= fBound) {
		pathLength = g;
		return true;
	}
	
	double fPossibleBound = 0;
    if (front2frontH){
		double h = env->HCost(currState, possibleMidState);
		fPossibleBound = std::max(fPossibleBound, g + h + possibleMidStateG);
		if (fgreater(fPossibleBound, fBound) || g > backwardBound){
			UpdateNextBound(fBound, fPossibleBound,backwardLoc);
			return false;
		}
	}
	if (isConsistent){
		double error = g - env->HCost(currState,originGoal);
		fPossibleBound = std::max(possibleMidStateG + otherH + error, fPossibleBound); 
		if (fgreater(fPossibleBound, fBound) || g > backwardBound){
			UpdateNextBound(fBound, fPossibleBound,backwardLoc);
			return false;
		}
	}
  
  double originalH = env->HCost(currState, originStart);
  fPossibleBound = std::max(g + originalH + otherError, fPossibleBound);
	if (fgreater(fPossibleBound, fBound) || g > backwardBound){
		UpdateNextBound(fBound, fPossibleBound,backwardLoc);
		return false;
	}
  
	std::vector<state> neighbors;
	env->GetSuccessors(currState, neighbors);
	nodesTouched += neighbors.size();
	nodesExpanded++;
	if(fPossibleBound == fBound){
		dMMLastIterExpansions++;
	}
	for (unsigned int x = 0; x < neighbors.size(); x++){
		uint64_t childID;
		double edgeCost = env->GCost(currState, neighbors[x]);
		if (neighbors[x] == parent || (detectDuplicates && backwardList.getElements().size()>1 && 
     backwardList.Lookup(env->GetStateHash(neighbors[x]), childID) != kNotFound && (backwardList.Lookup(childID).where == kClosedList || (backwardList.Lookup(childID).where == kOpenList && backwardList.Lookup(childID).g <= g + edgeCost)))) {
			continue;
		}
		if (DoIterationBackward(env, currState, neighbors[x], g + edgeCost, midState, possibleMidState, possibleMidStateG, otherH, otherError)) {
			return true;
		}
	}
	return false;
}

template <class state, class action, bool verbose>
void IDMM<state, action, verbose>::UpdateNextBound(double currBound, double fCost, int loc)
{
	if (!fgreater(nextBound[loc], currBound))
	{
		nextBound[loc] = fCost;
	}
	else if (fgreater(fCost, currBound) && fless(fCost, nextBound[loc]))
	{
		nextBound[loc] = fCost;
	}
}

#endif