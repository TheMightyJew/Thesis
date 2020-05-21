/*
 *  FullMBBDS.h
 *
 *  Created by Steven Danishevski on DEC 19.
 *
 */

#ifndef FullMBBDS_H
#define FullMBBDS_H

#include <iostream>
#include<algorithm> 
#include <unordered_set>
#include "SearchEnvironment.h"
#include "PancakeHasher.h"
#include <math.h>
#include "MM.h"
#include "IDMM.h"

template <class state, class action, class environment, class BloomFilter, const int pancakes_num, bool verbose = true>
class FullMBBDS {
public:
	FullMBBDS(unsigned long statesQuantityBound) {
		this->statesQuantityBound = statesQuantityBound;
	}
	virtual ~FullMBBDS() {}
	bool solve(SearchEnvironment<state, action>* env, state fromState, state toState, state &midState, std::vector<state> &thePath, int secondsLimit=600, bool threePhase=false);
	double getPathLength()	{ return pathLength; }
	uint64_t GetNodesExpanded() { return nodesExpanded; }
	bool isThreePhase() { return threePhase; }
private:
	unsigned long nodesExpanded, nodesTouched, statesQuantityBound;
	double pathLength;
	unsigned long memoryBound;
	bool threePhase=false;
};

template <class state, class action, class environment, class BloomFilter, const int pancakes_num, bool verbose>
bool FullMBBDS<state, action, environment, BloomFilter, pancakes_num, verbose>::solve(SearchEnvironment<state, action>* env,
	state fromState, state toState, state &midState, std::vector<state> &thePath, int secondsLimit, bool threePhase)
{
	this->threePhase = threePhase;
	bool solved = false;
	unsigned long currentNodesExapanded = 0;
	double lastBound = 0;
	if(threePhase){
		MM<state, action, environment> mm;
		solved = mm.GetPath(env, fromState, toState, env, env, thePath, secondsLimit, statesQuantityBound);
		currentNodesExapanded += mm.GetNodesExpanded();
		lastBound = mm.getLastBound();
	}
	if(solved){
		pathLength = env->GetPathLength(thePath);
		nodesExpanded = currentNodesExapanded;
		return true;
	}
	else{
		MBBDS<state, action, BloomFilter, false> mbbds(statesQuantityBound) ;
		PancakePuzzleState<pancakes_num> midState;
		solved = mbbds.GetMidState(env, fromState, toState, midState, secondsLimit, lastBound);
		currentNodesExapanded += mbbds.GetNodesExpanded();
		lastBound = mbbds.getLastBound();
		if(solved){
			pathLength = mbbds.getPathLength();
			nodesExpanded = currentNodesExapanded;
			return true;
		}
		else{
			IDMM<state, action, false> idmm;
			solved = idmm.GetMidState(env, fromState, toState, midState, secondsLimit, lastBound);
			currentNodesExapanded += idmm.GetNodesExpanded();
			if(solved){
				pathLength = idmm.getPathLength();
				nodesExpanded = currentNodesExapanded;
				return true;
			}
			else{
				return false;
			} 
		}
	}
}

#endif