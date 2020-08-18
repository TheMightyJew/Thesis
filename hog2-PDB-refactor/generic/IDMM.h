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
	IDMM(bool front2frontH=false, bool isConsistent = false) { this->front2frontH = front2frontH; this->isConsistent = isConsistent;}
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
	bool DoIterationBackward(SearchEnvironment<state, action>* env, state parent, state currState, double g, state& midState, state possibleMidState, double otherH = 0,double otherError = 0);
	void updateBoundsG(double minOpenG, double maxOpenG);
	void buildMatrix(std::vector<AStarOpenClosedData<state>> &openList, std::vector<std::vector<AStarOpenClosedData<state>>> &matrix);
  double nextBound;
	void UpdateNextBound(double currBound, double fCost);
	state originGoal;
	state originStart;
	unsigned long dMMLastIterExpansions = 0;
	unsigned long necessaryExpansions = 0;
	unsigned long nodesExpandedSoFar = 0;
	unsigned long previousIterationExpansions = 0;
	
	AStarOpenClosed<state, MMCompare<state>, AStarOpenClosedData<state>> forwardList;
	AStarOpenClosed<state, MMCompare<state>, AStarOpenClosedData<state>> backwardList;
	std::vector<AStarOpenClosedData<state>> forwardOpenList;
	std::vector<AStarOpenClosedData<state>> backwardOpenList;
	std::vector<std::vector<AStarOpenClosedData<state>>> forwardMatrix;
	std::vector<std::vector<AStarOpenClosedData<state>>> backwardMatrix;
	bool readyOpenLists;
	bool front2frontH;
	bool firstBounds;
  bool isConsistent;
  double minForwardError = 0;
  double minBackwardError = 0;
  //double prevMinForwardError = 0;
  //double prevMinBackwardError = 0;

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
  nextBound = 0;
	double minF = std::numeric_limits<double>::max();
  double minFforward = std::numeric_limits<double>::max();
  double minFbackward = std::numeric_limits<double>::max();
	double maxOpenG = 0;
	double minOpenG = std::numeric_limits<double>::max();
  minForwardError = DBL_MAX;
  minBackwardError = DBL_MAX;
  
	for (int x = 0; x < forwardList.OpenSize(); x++){
		AStarOpenClosedData<state> openState = forwardList.getElements()[forwardList.GetOpenItem(x)];
		
		minOpenG = std::min(minOpenG, openState.g);
		maxOpenG = std::max(maxOpenG , openState.g);
    minFforward = std::min(minFforward,openState.g+openState.h);
    minForwardError = std::min(minForwardError, openState.g - env->HCost(originStart, openState.data));
		forwardOpenList.push_back(openState);
	}
	//buildMatrix(forwardOpenList, forwardMatrix);
	//minF = forwardOpenList.front().g + forwardOpenList.front().h;
	//if (!isConsistent || forwardError == DBL_MAX){
  //  forwardError = 0;
  //}
	for (int x = 0; x < backwardList.OpenSize(); x++){
		AStarOpenClosedData<state> openState = backwardList.getElements()[backwardList.GetOpenItem(x)];
    minFbackward = std::min(minFbackward,openState.g+openState.h);
    minBackwardError = std::min(minBackwardError, openState.g - env->HCost(openState.data,originGoal));
		backwardOpenList.push_back(openState);
	}
  if (!isConsistent){
    minBackwardError = minForwardError = 0;
  }
	//buildMatrix(backwardOpenList, backwardMatrix);
	
	//minF = std::min(forwardOpenList.front().g + forwardOpenList.front().h, backwardOpenList.front().g + backwardOpenList.front().h);
	//if (!isConsistent || backwardError == DBL_MAX){
  //  backwardError = 0;
  //}
  minF = std::max(minFbackward,minFforward);
  	sort( forwardOpenList.begin( ), forwardOpenList.end( ), [ ]( const AStarOpenClosedData<state>& lhs, const AStarOpenClosedData<state>& rhs )
	{
	   return lhs.h < rhs.h;
	});
  sort( backwardOpenList.begin( ), backwardOpenList.end( ), [ ]( const AStarOpenClosedData<state>& lhs, const AStarOpenClosedData<state>& rhs )
	{
	   return lhs.h < rhs.h;
	});
	double initialHeuristic = env->HCost(fromState, toState);
	startingFBound = std::max(startingFBound, std::max(initialHeuristic, minF));
	forwardBound = minOpenG;
	backwardBound = startingFBound - forwardBound;
	firstBounds = true;
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
    for (AStarOpenClosedData<state> openState: forwardOpenList){
      auto currentTime = std::chrono::steady_clock::now();
				std::chrono::duration<double> elapsed_seconds = currentTime-startTime;
				if(elapsed_seconds.count() >= secondsLimit){
					return false;
				}
				if(openState.g + openState.h > forwardBound+backwardBound || openState.g > forwardBound)
					continue;
		//for(std::vector<AStarOpenClosedData<state>> fVector:forwardMatrix){
		//	if(fVector.front().g+fVector.front().h > forwardBound+backwardBound){
		//		break;
		//	}
		//	for (AStarOpenClosedData<state> openState: fVector){
		//		auto currentTime = std::chrono::steady_clock::now();
		//		std::chrono::duration<double> elapsed_seconds = currentTime-startTime;
		//		if(elapsed_seconds.count() >= secondsLimit){
		//			return false;
		//		}
		//		if(openState.g > forwardBound)
		//			break;
				if(firstBounds || openState.g == forwardBound){
					solved = DoIterationForward(env, openState.data, openState.data, openState.g, midState);
				}
				if(solved){
					break;
				}
			}
			//if(solved){
			//	break;
			//}
		//}
		if(verbose){
			printf("Nodes expanded: %d(%d)\n", nodesExpanded-nodesExpandedSoFar, nodesExpanded);
		}
		if (solved) {
			dMMExpansions = previousIterationExpansions + dMMLastIterExpansions;
			pathLength = std::min(pathLength, backwardBound + forwardBound);

			necessaryExpansions += nodesExpandedSoFar;
			for (AStarOpenClosedData<state> forwardState : forwardList.getElements()){
				if(forwardState.where == kClosedList && forwardState.g+forwardState.h<pathLength){
					necessaryExpansions++;
				}
			}
			for (AStarOpenClosedData<state> backwardState : backwardList.getElements()){
				if(backwardState.where == kClosedList && backwardState.g+backwardState.h<pathLength){
					necessaryExpansions++;
				}
			}
			return true;
		}
		else{
			updateBoundsG(minOpenG, maxOpenG);			
		}
	}
	return false;
}

template <class state, class action, bool verbose>
void IDMM<state, action, verbose>::updateBoundsG(double minOpenG, double maxOpenG){
	if(maxOpenG > forwardBound && backwardBound > 0){
		forwardBound++;
		backwardBound--;
		firstBounds = false;
	}
	else{
		double fullBound = forwardBound + backwardBound + 1;
		forwardBound = std::max(minOpenG, ceil(fullBound/2));
		backwardBound = fullBound - forwardBound;
		firstBounds = true;
		previousIterationExpansions = nodesExpanded-nodesExpandedSoFar;
		nodesExpandedSoFar = nodesExpanded;
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
  nextBound = 0;
	nodesExpanded = nodesTouched = 0;
	originStart = fromState;
	originGoal = toState;
	double initialHeuristic = env->HCost(fromState, toState);
  UpdateNextBound(0,env->HCost(fromState, toState));
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
      forwardBound = ceil(nextBound / 2);
      backwardBound = nextBound - forwardBound;
      /*
			if (forwardBound > backwardBound) {
				backwardBound = forwardBound;
			}
			else{
				forwardBound++;
			}
      */
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
	if (g > forwardBound  || fgreater(g + h + minBackwardError, backwardBound + forwardBound))//not sure about that
	{
    UpdateNextBound(backwardBound + forwardBound,g + h + minBackwardError);
		return false;
	}
	else if (g == forwardBound) {
    double error = 0;
    if (isConsistent){
      error = g - env->HCost(currState,originStart);
    }
		if(readyOpenLists){	
			//for(std::vector<AStarOpenClosedData<state>> fVector:backwardMatrix){
			//	if(fVector.front().g+fVector.front().h > forwardBound+backwardBound){
			//		break;
			//	}
			//	for (AStarOpenClosedData<state> openState: fVector){
			//		if(openState.g > backwardBound)
			//			break;
          for (AStarOpenClosedData<state> openState: backwardOpenList){
            if(openState.g+openState.h > forwardBound+backwardBound || openState.g > backwardBound){
              UpdateNextBound(backwardBound + forwardBound,openState.g+openState.h);
              continue;
            }
            if (DoIterationBackward(env, openState.data, openState.data, openState.g, midState, currState,h,error)) {
              pathLength += g;
              midState = currState;
              return true;
            }
          }
				//}
			//}
			return false;
		}
		else if (DoIterationBackward(env, originGoal, originGoal, 0, midState, currState,h,error)) {
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
		double edgeCost = env->GCost(currState, neighbors[x]);
		if (neighbors[x] == parent || (readyOpenLists && forwardList.Lookup(env->GetStateHash(neighbors[x]), childID) != kNotFound && (forwardList.Lookup(childID).where == kClosedList || (forwardList.Lookup(childID).where == kOpenList && forwardList.Lookup(childID).g <= g + edgeCost)))) {
			continue;
		}
		/*if (neighbors[x] == parent || forwardList.Lookup(env->GetStateHash(neighbors[x]), childID) == kClosedList) {
			continue;
		}*/
		if (DoIterationForward(env, currState, neighbors[x], g + edgeCost, midState)) {
			return true;
		}
	}
	return false;
}

template <class state, class action, bool verbose>
bool IDMM<state, action, verbose>::DoIterationBackward(SearchEnvironment<state, action>* env,
	state parent, state currState, double g, state& midState, state possibleMidState, double otherH,double otherError)
{	
	if (currState == possibleMidState) {
		pathLength = g;
		return true;
	}
    double fbound = 0;
  	if (g>=backwardBound){
    UpdateNextBound(forwardBound+backwardBound, fbound);
		return false;
	} 
    if (front2frontH){
    double h = env->HCost(currState, possibleMidState);
    fbound = std::max(fbound, g + h + forwardBound);
    if (fgreater(fbound, forwardBound+backwardBound)){
      UpdateNextBound(forwardBound+backwardBound, fbound);
      return false;
    }
  }
  if (isConsistent){
    double error = g - env->HCost(currState,originGoal);
    fbound = std::max(forwardBound + otherH + error,fbound); 
    if (fgreater(fbound, forwardBound+backwardBound)){
      UpdateNextBound(forwardBound+backwardBound, fbound);
      return false;
    }
  }

  double originalH = env->HCost(currState, originStart);
  fbound = std::max(g + originalH + otherError,fbound);
  if (fgreater(fbound, forwardBound+backwardBound) || g>=backwardBound){
    UpdateNextBound(forwardBound+backwardBound, fbound);
    return false;
  }
  
	std::vector<state> neighbors;
	env->GetSuccessors(currState, neighbors);
	nodesTouched += neighbors.size();
	nodesExpanded++;
	if(fbound == forwardBound+backwardBound){
		dMMLastIterExpansions++;
	}
	for (unsigned int x = 0; x < neighbors.size(); x++)
	{
		uint64_t childID;
		double edgeCost = env->GCost(currState, neighbors[x]);
		if (neighbors[x] == parent || (readyOpenLists && (backwardList.Lookup(env->GetStateHash(neighbors[x]), childID) == kClosedList || (backwardList.Lookup(env->GetStateHash(neighbors[x]), childID) == kOpenList && backwardList.Lookup(childID).g <= g + edgeCost)))) {
			continue;
		}
		/*if (neighbors[x] == parent || backwardList.Lookup(env->GetStateHash(neighbors[x]), childID) == kClosedList) {
			continue;
		}*/
		if (DoIterationBackward(env, currState, neighbors[x], g + edgeCost, midState, possibleMidState,otherH,otherError)) {
			return true;
		}
	}
	return false;
}

template <class state, class action, bool verbose>
void IDMM<state, action, verbose>::UpdateNextBound(double currBound, double fCost)
{
	if (!fgreater(nextBound, currBound))
	{
		nextBound = fCost;
		//printf("Updating next bound to %f\n", nextBound);
	}
	else if (fgreater(fCost, currBound) && fless(fCost, nextBound))
	{
		nextBound = fCost;
		//printf("Updating next bound to %f\n", nextBound);
	}
}

#endif