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
#include "AStarOpenClosed.h"

//#define DO_LOGGING

typedef __gnu_cxx::hash_map<uint64_t, double> NodeHashTable;

template <class state, class action, bool verbose = true>
class IDAStar {
public:
	IDAStar() { useHashTable = usePathMax = false; storedHeuristic = false;}
	virtual ~IDAStar() {}
	bool GetPath(SearchEnvironment<state, action> *env, state from, state to,
							 std::vector<state> &thePath, int secondsLimit=600, bool readyOpenList=false, std::vector<AStarOpenClosedDataWithF<state>> openList=std::vector<AStarOpenClosedDataWithF<state>>());
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
					   double maxH, double currStateH=0);
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
	double nextBound;
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

#ifdef DO_LOGGING
public:
	std::function<void (state, int)> func;
#endif
};

template <class state, class action, bool verbose>
bool IDAStar<state, action, verbose>::GetPath(SearchEnvironment<state, action> *env,
									 state from, state to,
									 std::vector<state> &thePath, int secondsLimit, bool readyOpenList, std::vector<AStarOpenClosedDataWithF<state>> openList)
{
	if(verbose){
		printf("\t\tStarting to solve with IDAStar\n");
	}
	if (!storedHeuristic)
		heuristic = env;
	nextBound = 0;
	nodesExpanded = nodesTouched = dAstarLastIterExpansions = 0;
	thePath.resize(0);
	if(readyOpenList){
		double minF = 0;
		double maxF = 0;
		for(AStarOpenClosedDataWithF<state> openState : openList){
			if(minF == 0 || minF > openState.f){
					minF = openState.f;
			}
			if(maxF == 0 || maxF < openState.f){
					maxF = openState.f;
			}
		}
		UpdateNextBound(0, minF);
		
		for(AStarOpenClosedDataWithF<state> openState : openList){
			if(maxF > openState.f){
				necessaryExpansions++;
			}
		}
	}
	else{
		UpdateNextBound(0, heuristic->HCost(from, to));
	}
	goal = to;
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
		double res;
		if(readyOpenList){
			double currentBound = nextBound;
			for(AStarOpenClosedDataWithF<state> openState : openList){
				res = DoIteration(env, openState.data, openState.data, thePath, currentBound, openState.g, openState.h);
				if(res == 0 && solved){
					solLength = env->GetPathLength(thePath) + openState.g;
					break;
				}
			}
		}
		else{
			res = DoIteration(env, from, from, thePath, nextBound, 0, 0);
		}
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
										   double maxH, double currStateH)
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
	if (env->GoalTest(currState, goal)){
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

	for (unsigned int x = 0; x < neighbors.size(); x++)
	{
		if (neighbors[x] == parent)
			continue;
		thePath.push_back(neighbors[x]);
		double edgeCost = env->GCost(currState, neighbors[x]);
		double childH = DoIteration(env, currState, neighbors[x], thePath, bound,
																g+edgeCost, maxH - edgeCost);
		if (env->GoalTest(thePath.back(), goal) && g+edgeCost<=bound){
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

