/*
 *  IDAStar.h
 *  hog2
 *
 *  Created by Nathan Sturtevant on 5/22/07.
 *  Copyright 2007 Nathan Sturtevant, University of Alberta. All rights reserved.
 *
 */

#ifndef IDASTAR_H
#define IDASTAR_H

#include <iostream>
#include <functional>
#include "SearchEnvironment.h"
#include <ext/hash_map>
#include "FPUtil.h"
#include "vectorCache.h"
#include "TemplateAStar.h"
#include "AStarOpenClosed.h"
#include <limits>
#include <algorithm>

//#define DO_LOGGING

typedef __gnu_cxx::hash_map<uint64_t, double> NodeHashTable;

template <class state, class action, bool verbose = true>
class IDAStar {
public:
	IDAStar() { useHashTable = usePathMax = false; storedHeuristic = false;}
	virtual ~IDAStar() {}
	bool GetPath(SearchEnvironment<state, action> *env, state from, state to,
							 std::vector<state> &thePath, int secondsLimit=600);
	bool ASpIDA(SearchEnvironment<state, action> *env, state from, state to, std::vector<state> &thePath, AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> statesList, int secondsLimit=600);
	bool BAI(SearchEnvironment<state, action> *env, state from, state to, std::vector<state> &thePath, AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> statesList, int secondsLimit=600, bool isConsistent = false);
  bool ASpIDArev(SearchEnvironment<state, action> *env, state from, state to, std::vector<state> &thePath, AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> statesList, int secondsLimit=600, bool isConsistent = false);
	bool GetPath(SearchEnvironment<state, action> *env, state from, state to,
				 std::vector<action> &thePath);

	uint64_t GetNodesExpanded() { return nodesExpanded; }
	uint64_t GetNecessaryExpansions() { return necessaryExpansions; }
	uint64_t GetNodesTouched() { return nodesTouched; }
	void ResetNodeCount() { nodesExpanded = nodesTouched = 0; }
	void SetUseBDPathMax(bool val) { usePathMax = val; }
	void SetHeuristic(Heuristic<state> *heur) { heuristic = heur; if (heur != 0) storedHeuristic = true;}
	unsigned long getDAstarExpansions() { return dAstarExpansions; }
	double getSolLength() { return solLength; }
private:
	unsigned long long nodesExpanded, nodesTouched;
	
	double DoIteration(SearchEnvironment<state, action> *env,
					   state parent, state currState,
					   std::vector<state> &thePath, double bound, double g,
					   double maxH, double currStateH=0, bool updateH=false, AStarOpenClosedDataWithF<state> *stateObject=NULL);
	double DoIteration(SearchEnvironment<state, action> *env,
					   action forbiddenAction, state &currState,
					   std::vector<action> &thePath, double bound, double g,
					   double maxH, double parentH);
	void PrintGHistogram()
	{
//		uint64_t early = 0, late = 0;
//		for (int x = 0; x < gCostHistogram.size(); x++)
//		{
//			printf("%d\t%llu\n", x, gCostHistogram[x]);
//			if (x*2 > gCostHistogram.size()-1)
//				late += gCostHistogram[x];
//			else
//				early += gCostHistogram[x];
//		}
//		if (late < early)
//			printf("Strong heuristic - Expect MM > A*\n");
//		else
//			printf("Weak heuristic - Expect MM >= MM0.\n");
	}
	void UpdateNextBound(double currBound, double fCost);
	state goal;
  state start;
	double nextBound, currIterNextBound;
	//NodeHashTable nodeTable;
	bool usePathMax;
	bool useHashTable;
	vectorCache<action> actCache;
	bool storedHeuristic;
	bool solved = false;
	double solLength;
	Heuristic<state> *heuristic;
	std::vector<uint64_t> gCostHistogram;
	unsigned long dAstarExpansions = 0;
	unsigned long dAstarLastIterExpansions = 0;
	unsigned long necessaryExpansions = 0;
	AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> statesList;
	std::vector<AStarOpenClosedDataWithF<state>> openList;
	//std::vector<AStarOpenClosedDataWithF<state>> perimeterList;
	double perimeterG;
	bool readyStatesList;
	bool reverseG = false;
	bool reverseF = false;
	//double minError;
  bool isConsistent = false;
  double minPerimeterF;

#ifdef DO_LOGGING
public:
	std::function<void (state, int)> func;
#endif
};

template <class state, class action, bool verbose>
bool IDAStar<state, action, verbose>::ASpIDA(SearchEnvironment<state, action> *env, state from, state to,
							 std::vector<state> &thePath, AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> statesList, int secondsLimit)
{
	//reverseG = false;
	heuristic = env;
	nextBound = 0;
	nodesExpanded = nodesTouched = dAstarLastIterExpansions = 0;
	thePath.resize(0);
	this->readyStatesList = true;
	this->statesList = statesList;
	double minF = std::numeric_limits<double>::max();
	double maxF = 0;
	for (int x = 0; x < statesList.OpenSize(); x++){
		AStarOpenClosedDataWithF<state> openState = statesList.getElements()[statesList.GetOpenItem(x)];
		openList.push_back(openState);
		minF = std::min(minF, openState.f);
		maxF = std::max(maxF, openState.f);
	}
	UpdateNextBound(0, std::max(minF, heuristic->HCost(from, to)));
	sort( openList.begin( ), openList.end( ), [ ]( const AStarOpenClosedDataWithF<state>& lhs, const AStarOpenClosedDataWithF<state>& rhs )
	{
	   return lhs.h < rhs.h;
	});
	goal = to;
	start = from;
	thePath.push_back(from);
	unsigned long nodesExpandedSoFar = 0;
	unsigned long previousIterationExpansions = 0;
	auto startTime = std::chrono::steady_clock::now();
	while (true) //thePath.size() == 0)
	{
		gCostHistogram.clear();
		gCostHistogram.resize(nextBound+1);
		dAstarLastIterExpansions = 0;
		auto currentTime = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = currentTime-startTime;
		if(elapsed_seconds.count() >= secondsLimit){
			return false;
		}
		if (verbose)
			printf("\t\tStarting iteration with bound %1.1f: ", nextBound, nodesExpanded);
		double res;
		double currentBound = nextBound;
		for(AStarOpenClosedDataWithF<state> openState:openList){
			thePath.resize(0);
			thePath.push_back(openState.data);
			currIterNextBound = currentBound;
			res = DoIteration(env, openState.data, openState.data, thePath, currentBound, openState.g, 0, openState.h, true, &openState);
			if(res == 0 && solved){
				solLength = env->GetPathLength(thePath) + openState.g;
				break;
			}
			else{	
				openState.h += currIterNextBound -(openState.g+openState.h);	
            }
		}
		if (verbose)
			printf("Nodes expanded: %llu(%llu)\n", nodesExpanded-nodesExpandedSoFar, nodesExpanded);
		if (res == 0 && solved){
			if (verbose){
				printf("\t\tStarting iteration with bound %f. ", nextBound, nodesExpanded);
				printf("Nodes expanded: %llu(%llu)\n", nodesExpanded-nodesExpandedSoFar, nodesExpanded);
			}
			necessaryExpansions += nodesExpandedSoFar;
			for (AStarOpenClosedDataWithF<state> astarState : statesList.getElements()){
				if(astarState.where == kClosedList && astarState.f<solLength){
					necessaryExpansions++;
				}
			}
			break;
		}
		previousIterationExpansions = nodesExpanded-nodesExpandedSoFar;
		nodesExpandedSoFar = nodesExpanded;
		PrintGHistogram();
	}
	dAstarExpansions = previousIterationExpansions + dAstarLastIterExpansions;
	PrintGHistogram();
	return true;
}

template <class state, class action, bool verbose>
bool IDAStar<state, action, verbose>::BAI(SearchEnvironment<state, action> *env, state from, state to,
							 std::vector<state> &thePath, AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> statesList, int secondsLimit, bool isConsistent)
{
  this->isConsistent = isConsistent;
	reverseF = true;
	heuristic = env;
	nextBound = 0;
	nodesExpanded = nodesTouched = dAstarLastIterExpansions = 0;
	thePath.resize(0);
  this->statesList = statesList;
	this->readyStatesList = false;
	//double minF = std::numeric_limits<double>::max();
	//for (AStarOpenClosedDataWithF<state> astarState : statesList.getElements()){
		//minF = std::min(minF, astarState.g + astarState.h);
		//maxF = std::max(maxF, astarState.f);
		//maxG = std::max(maxG, astarState.g);
		//if (isConsistent){
		//	minError = std::min(minError, astarState.g - heuristic->HCost(astarState.data, to));
		//}
	//}
  uint64_t key = statesList.Peek();
	minPerimeterF = statesList.Lookup(key).g+statesList.Lookup(key).h;

	UpdateNextBound(0, std::max(minPerimeterF, heuristic->HCost(from, to)));
	goal = to;
	start = from;
	thePath.push_back(from);
	unsigned long nodesExpandedSoFar = 0;
	unsigned long previousIterationExpansions = 0;
	auto startTime = std::chrono::steady_clock::now();
	while (true)
	{
		gCostHistogram.clear();
		gCostHistogram.resize(nextBound+1);
		dAstarLastIterExpansions = 0;
		auto currentTime = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = currentTime-startTime;
		if(elapsed_seconds.count() >= secondsLimit){
			return false;
		}
		if (verbose)
			printf("\t\tStarting iteration with bound %1.1f: ", nextBound, nodesExpanded);
		double res = DoIteration(env, from, from, thePath, nextBound, 0, 0);
		if (verbose)
			printf("Nodes expanded: %llu(%llu)\n", nodesExpanded-nodesExpandedSoFar, nodesExpanded);
		if (res == 0 && solved){
			if (verbose){
				printf("\t\tStarting iteration with bound %f. ", nextBound, nodesExpanded);
				printf("Nodes expanded: %llu(%llu)\n", nodesExpanded-nodesExpandedSoFar, nodesExpanded);
			}
			necessaryExpansions += nodesExpandedSoFar;
			for (AStarOpenClosedDataWithF<state> astarState : statesList.getElements()){
				if(astarState.where == kClosedList && astarState.f<solLength){
					necessaryExpansions++;
				}
			}
			break;
		}
		previousIterationExpansions = nodesExpanded-nodesExpandedSoFar;
		nodesExpandedSoFar = nodesExpanded;
		PrintGHistogram();
	}
	dAstarExpansions = previousIterationExpansions + dAstarLastIterExpansions;
	PrintGHistogram();
	return true;
}

template <class state, class action, bool verbose>
bool IDAStar<state, action, verbose>::ASpIDArev(SearchEnvironment<state, action> *env, state from, state to,
							 std::vector<state> &thePath, AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> statesList, int secondsLimit, bool isConsistent)
{
	reverseG = true;
	heuristic = env;
	nextBound = 0;
	nodesExpanded = nodesTouched = dAstarLastIterExpansions = 0;
	thePath.resize(0);
  this->isConsistent = isConsistent;
	this->readyStatesList = false;
  this->statesList = statesList;
	//double maxF = 0;
	double minRealF = std::numeric_limits<double>::max();
	double maxG = 0;
	for (AStarOpenClosedDataWithF<state> astarState : statesList.getElements()){
		minRealF = std::min(minRealF, astarState.g + heuristic->HCost(astarState.data, from));
		//maxF = std::max(maxF, astarState.f);
		//maxG = std::max(maxG, astarState.g);
	}
  uint64_t key = statesList.Peek();
	perimeterG = std::max(0.0, statesList.Lookup(key).f - 1); // this needs to be fixed for non-integer domains
  /*
	for (AStarOpenClosedDataWithF<state> astarState : statesList.getElements()){
		//changed
		if(astarState.g < maxG){
			perimeterList.push_back(astarState);
		}
	}
  */
	UpdateNextBound(0, std::max(minRealF, heuristic->HCost(from, to)));
	goal = to;
	start = from;
	thePath.push_back(from);
	unsigned long nodesExpandedSoFar = 0;
	unsigned long previousIterationExpansions = 0;
	auto startTime = std::chrono::steady_clock::now();
	while (true)
	{
		gCostHistogram.clear();
		gCostHistogram.resize(nextBound+1);
		dAstarLastIterExpansions = 0;
		auto currentTime = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = currentTime-startTime;
		if(elapsed_seconds.count() >= secondsLimit){
			return false;
		}
		if (verbose)
			printf("\t\tStarting iteration with bound %1.1f: ", nextBound, nodesExpanded);
		double res = DoIteration(env, from, from, thePath, nextBound, 0, 0);
		if (verbose)
			printf("Nodes expanded: %llu(%llu)\n", nodesExpanded-nodesExpandedSoFar, nodesExpanded);
		if (res == 0 && solved){
			if (verbose){
				printf("\t\tStarting iteration with bound %f. ", nextBound, nodesExpanded);
				printf("Nodes expanded: %llu(%llu)\n", nodesExpanded-nodesExpandedSoFar, nodesExpanded);
			}
			necessaryExpansions += nodesExpandedSoFar;
			for (AStarOpenClosedDataWithF<state> astarState : statesList.getElements()){
				if(astarState.where == kClosedList && astarState.f<solLength){
					necessaryExpansions++;
				}
			}
			break;
		}
		previousIterationExpansions = nodesExpanded-nodesExpandedSoFar;
		nodesExpandedSoFar = nodesExpanded;
		PrintGHistogram();
	}
	dAstarExpansions = previousIterationExpansions + dAstarLastIterExpansions;
	PrintGHistogram();
	return true;
}

template <class state, class action, bool verbose>
bool IDAStar<state, action, verbose>::GetPath(SearchEnvironment<state, action> *env,
									 state from, state to,
									 std::vector<state> &thePath, int secondsLimit)
{
	if(verbose){
		printf("\t\tStarting to solve with IDAStar\n");
	}
	if (!storedHeuristic)
		heuristic = env;
	nextBound = 0;
	nodesExpanded = nodesTouched = dAstarLastIterExpansions = 0;
	thePath.resize(0);
	this->readyStatesList = false;
	UpdateNextBound(0, heuristic->HCost(from, to));
	goal = to;
	start = from;
	thePath.push_back(from);
	unsigned long nodesExpandedSoFar = 0;
	unsigned long previousIterationExpansions = 0;
	auto startTime = std::chrono::steady_clock::now();
	while (true) //thePath.size() == 0)
	{
		dAstarLastIterExpansions = 0;
		auto currentTime = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = currentTime-startTime;
		if(elapsed_seconds.count() >= secondsLimit){
			return false;
		}
		//nodeTable.clear();
		gCostHistogram.clear();
		gCostHistogram.resize(nextBound+1);
		if (verbose)
			printf("\t\tStarting iteration with bound %1.1f: ", nextBound, nodesExpanded);
		double res = DoIteration(env, from, from, thePath, nextBound, 0, 0);
		if (verbose)
			printf("Nodes expanded: %llu(%llu)\n", nodesExpanded-nodesExpandedSoFar, nodesExpanded);
		if (res == 0 && solved){
			if (verbose){
				printf("\t\tStarting iteration with bound %f. ", nextBound, nodesExpanded);
				printf("Nodes expanded: %llu(%llu)\n", nodesExpanded-nodesExpandedSoFar, nodesExpanded);
			}
			necessaryExpansions += nodesExpandedSoFar;
			break;
		}
		previousIterationExpansions = nodesExpanded-nodesExpandedSoFar;
		nodesExpandedSoFar = nodesExpanded;
		PrintGHistogram();
	}
	dAstarExpansions = previousIterationExpansions + dAstarLastIterExpansions;
	PrintGHistogram();
	return true;
}

template <class state, class action, bool verbose>
bool IDAStar<state, action, verbose>::GetPath(SearchEnvironment<state, action> *env,
									 state from, state to,
									 std::vector<action> &thePath)
{
	if (!storedHeuristic)
		heuristic = env;
	nextBound = 0;
	nodesExpanded = nodesTouched = 0;
	thePath.resize(0);

	if (env->GoalTest(from, to))
		return true;

	double rootH = heuristic->HCost(from, to);
	UpdateNextBound(0, rootH);
	goal = to;
	start = from;
	std::vector<action> act;
	env->GetActions(from, act);
	while (thePath.size() == 0)
	{
		//nodeTable.clear();
		gCostHistogram.clear();
		gCostHistogram.resize(nextBound+1);
		if (verbose)
			printf("Starting iteration with bound %f; %llu expanded, %llu generated\n", nextBound, nodesExpanded, nodesTouched);
		fflush(stdout);
		DoIteration(env, act[0], from, thePath, nextBound, 0, 0, rootH);
		PrintGHistogram();
	}
	return true;
}

template <class state, class action, bool verbose>
double IDAStar<state, action, verbose>::DoIteration(SearchEnvironment<state, action> *env,
										   state parent, state currState,
										   std::vector<state> &thePath, double bound, double g,
										   double maxH, double currStateH, bool updateH, AStarOpenClosedDataWithF<state> *stateObject)
{
	double h = std::max(heuristic->HCost(currState, goal), currStateH);
	// path max
	if (usePathMax && fless(h, maxH))
		h = maxH;
	if (fgreater(g+h, bound))
	{
		UpdateNextBound(bound, g+h);
		//printf("Stopping at (%d, %d). g=%f h=%f\n", currState>>16, currState&0xFFFF, g, h);
		return h;
	}

	else if(reverseG && perimeterG+g >= bound){

    uint64_t ID;
    if (statesList.Lookup(env->GetStateHash(currState), ID) != kNotFound){
      if (statesList.Lookup(ID).g + g <= bound){
				solLength = g + statesList.Lookup(ID).g;
				solved = true;
				return 0;
      }
      else{
        UpdateNextBound(bound, statesList.Lookup(ID).g + g);
        return h;
      }
		}
    else if (fgreater(perimeterG+g,bound)){
      UpdateNextBound(bound, g+perimeterG);
      return h;
    }
	}
  else if(reverseF){
    uint64_t ID;
    if (statesList.Lookup(env->GetStateHash(currState), ID) != kNotFound && statesList.Lookup(ID).g + g <= bound){
				solLength = g + statesList.Lookup(ID).g;
				solved = true;
				return 0;
		}
    else if (isConsistent){
      double error = g - heuristic->HCost(currState, start);
      if (fgreater(minPerimeterF + error ,bound)){
        UpdateNextBound(bound, minPerimeterF + error);
        return h; 
      }      
    }
	}
  
  else if (env->GoalTest(currState, goal)){
		solved = true;
		return 0;
	}
  		
	std::vector<state> neighbors;
	env->GetSuccessors(currState, neighbors);
	nodesTouched += neighbors.size();
	nodesExpanded++;
	if(g+h == bound){
		dAstarLastIterExpansions++;
	}
	gCostHistogram[g]++;

	double minNeighborF = std::numeric_limits<double>::max();
	for (unsigned int x = 0; x < neighbors.size(); x++){
		uint64_t childID;
		if (neighbors[x] == parent || (readyStatesList && statesList.Lookup(env->GetStateHash(neighbors[x]), childID) == kClosedList)) {
			continue;
		}
		thePath.push_back(neighbors[x]);
		double edgeCost = env->GCost(currState, neighbors[x]);
		double childH = DoIteration(env, currState, neighbors[x], thePath, bound,
																g+edgeCost, maxH - edgeCost);
		/*if (updateH){
		  minNeighborF = std::min(minNeighborF, g+edgeCost+heuristic->HCost(neighbors[x], goal)+minError);	
		}*/
		if(solved){
			return 0;
		}
		else if (env->GoalTest(thePath.back(), goal) && g+edgeCost<=bound){
			solved = true;
			return 0;
		}
		thePath.pop_back();
		// pathmax
		if (usePathMax && fgreater(childH-edgeCost, h))
		{
//			nodeTable[currState] = g;//+h
			h = childH-edgeCost;
			if (fgreater(g+h, bound))
			{
				UpdateNextBound(bound, g+h);
				return h;
			}
		}
	}
	/*if(updateH){
		stateObject->h = stateObject->h + (minNeighborF-g);
	}*/
	return h;
}

template <class state, class action, bool verbose>
double IDAStar<state, action, verbose>::DoIteration(SearchEnvironment<state, action> *env,
										   action forbiddenAction, state &currState,
										   std::vector<action> &thePath, double bound, double g,
										   double maxH, double parentH)
{
	double h = heuristic->HCost(currState, goal);//, parentH); // TODO: restore code that uses parent h-cost
	parentH = h;
	// path max
	if (usePathMax && fless(h, maxH))
		h = maxH;
	if (fgreater(g+h, bound))
	{
		UpdateNextBound(bound, g+h);
		//printf("Stopping at (%d, %d). g=%f h=%f\n", currState>>16, currState&0xFFFF, g, h);
		return h;
	}
	// must do this after we check the f-cost bound
	if (env->GoalTest(currState, goal))
		return -1; // found goal
	
	std::vector<action> &actions = *actCache.getItem();
	env->GetActions(currState, actions);
	nodesTouched += actions.size();
	nodesExpanded++;
	gCostHistogram[g]++;
	int depth = (int)thePath.size();
#ifdef t
	func(currState, depth);
#endif
	
	for (unsigned int x = 0; x < actions.size(); x++)
	{
		if ((depth != 0) && (actions[x] == forbiddenAction))
			continue;

		thePath.push_back(actions[x]);
		double edgeCost = env->GCost(currState, actions[x]);
		env->ApplyAction(currState, actions[x]);
		action a = actions[x];
		env->InvertAction(a);

		double childH = DoIteration(env, a, currState, thePath, bound,
									g+edgeCost, maxH - edgeCost, parentH);
		env->UndoAction(currState, actions[x]);
		if (fequal(childH, -1)) // found goal
		{
			actCache.returnItem(&actions);
			return -1;
		}

		thePath.pop_back();

		// pathmax
		if (usePathMax && fgreater(childH-edgeCost, h))
		{
			//			nodeTable[currState] = g;//+h
			h = childH-edgeCost;
			if (fgreater(g+h, bound))
			{
				UpdateNextBound(bound, g+h);
				actCache.returnItem(&actions);
				return h;
			}
		}
	}
	actCache.returnItem(&actions);
	return h;
}


template <class state, class action, bool verbose>
void IDAStar<state, action, verbose>::UpdateNextBound(double currBound, double fCost)
{
	if (!fgreater(nextBound, currBound))
	{
		nextBound = currIterNextBound = fCost;
		//printf("Updating next bound to %f\n", nextBound);
	}
	else if (fgreater(fCost, currBound) && fless(fCost, nextBound))
	{
		nextBound = currIterNextBound = fCost;
		//printf("Updating next bound to %f\n", nextBound);
	}
}


#endif

