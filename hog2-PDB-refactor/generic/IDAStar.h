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
	bool ASpIDA(SearchEnvironment<state, action> *env, state from, state to, std::vector<state> &thePath, AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> statesList, int secondsLimit=600, bool isDuplicateDetection = true);
	bool BAI(SearchEnvironment<state, action> *env, state from, state to, std::vector<state> &thePath, AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> statesList, int secondsLimit=600, bool isConsistent = false);
  bool ASpIDArev(SearchEnvironment<state, action> *env, state from, state to, std::vector<state> &thePath, AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> statesList, int secondsLimit=600, bool isConsistent = false, bool isComputeMaxH = false);
	bool GetPath(SearchEnvironment<state, action> *env, state from, state to,
				 std::vector<action> &thePath);
         
  bool SFBDS(SearchEnvironment<state, action> *env, state from, state to,
							 std::vector<state> &thePath, int secondsLimit, int lookahead,bool isAlternating, bool storeLookup);

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
  bool DoSFBDSIteration(SearchEnvironment<state, action> *env,
             state parentStart, state currStart, state parentGoal, state currGoal,
             std::vector<state> &startPath, std::vector<state> &goalPath, double bound, double g_from, double g_to, int k,uint64_t b,double h, bool isAlternating,bool storeLookup, bool isForward);
  double computeKlook(SearchEnvironment<state, action> *env,state parent, state curr, state opposite,double bound,double g_curr,double g_opposite,double prevf,int k,std::vector<state> &thePath,bool &isSolved, uint64_t b,bool storeLookup, std::vector<state> &stateStorage, std::vector<double> &hStorage);
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
	std::vector<uint64_t> heuristicList;
	AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> statesList;
	std::vector<AStarOpenClosedDataWithF<state>> openList;
	//std::vector<AStarOpenClosedDataWithF<state>> perimeterList;
	double perimeterG;
	bool isDuplicateDetection;
	bool reverseG = false;
	bool reverseF = false;
	//double minError;
  bool isConsistent = false;
  double minPerimeterF;
  bool isComputeMaxH = false;

#ifdef DO_LOGGING
public:
	std::function<void (state, int)> func;
#endif
};

template <class state, class action, bool verbose>
bool IDAStar<state, action, verbose>::ASpIDA(SearchEnvironment<state, action> *env, state from, state to,
							 std::vector<state> &thePath, AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> statesList, int secondsLimit, bool isDuplicateDetection)
{
	//reverseG = false;
	heuristic = env;
	nextBound = 0;
	nodesExpanded = nodesTouched = dAstarLastIterExpansions = 0;
	thePath.resize(0);
	this->isDuplicateDetection = isDuplicateDetection;
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
			res = DoIteration(env, openState.data, openState.data, thePath, currentBound, openState.g, 0, openState.h);
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
	this->isDuplicateDetection = false;
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
							 std::vector<state> &thePath, AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> statesList, int secondsLimit, bool isConsistent, bool isComputeMaxH)
{
	reverseG = true;
	heuristic = env;
	nextBound = 0;
	nodesExpanded = nodesTouched = dAstarLastIterExpansions = 0;
	thePath.resize(0);
  this->isConsistent = isConsistent;
  this->isComputeMaxH = isComputeMaxH;
	this->isDuplicateDetection = false;
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
	double h = heuristic->HCost(from, to);
	UpdateNextBound(0, std::max(minRealF, h));
	goal = to;
	start = from;
	thePath.push_back(from);
	unsigned long nodesExpandedSoFar = 0;
	unsigned long previousIterationExpansions = 0;
	auto startTime = std::chrono::steady_clock::now();
	while (true)
	{
		double bound = nextBound;
		gCostHistogram.clear();
		gCostHistogram.resize(bound+1);
		dAstarLastIterExpansions = 0;
		auto currentTime = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = currentTime-startTime;
		if(elapsed_seconds.count() >= secondsLimit){
			return false;
		}
		if (verbose)
			printf("\t\tStarting iteration with bound %1.1f: ", bound, nodesExpanded);
		
		
		if(isComputeMaxH && isConsistent){
			heuristicList.clear();
			for (AStarOpenClosedDataWithF<state> astarState : statesList.getElements()){
				if(astarState.g == perimeterG){
					double calculatedF = heuristic->HCost(from, astarState.data) + astarState.g;
					if(bound >= calculatedF){
						uint64_t childID;
						statesList.Lookup(env->GetStateHash(astarState.data), childID);
						heuristicList.push_back(childID);
					}
					else{
						UpdateNextBound(bound, calculatedF);
					}
				}
			}	
		}
		double res;
		if(isComputeMaxH && isConsistent && heuristicList.size() == 0)
			res = 1;
		else
			res = DoIteration(env, from, from, thePath, bound, 0, 0);
		if (verbose)
			printf("Nodes expanded: %llu(%llu)\n", nodesExpanded-nodesExpandedSoFar, nodesExpanded);
		if (res == 0 && solved){
			if (verbose){
				printf("\t\tStarting iteration with bound %f. ", bound, nodesExpanded);
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
	this->isDuplicateDetection = false;
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
										   double maxH, double currStateH)
{
	std::vector<uint64_t> ignoreList;
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
	else if (reverseG && isComputeMaxH){
		double F2F_h = DBL_MAX;
    if (isConsistent){
      for (int i = 0; i < heuristicList.size(); i++){
        uint64_t childID = heuristicList[i];
        AStarOpenClosedDataWithF<state> perimeterState = statesList.Lookup(childID);
        if(perimeterState.g == perimeterG){ //perimeter frontier should be fixed and this should be updated to <= instead 
          double calculatedF = g + heuristic->HCost(currState, perimeterState.data) + perimeterState.g;
          F2F_h = std::min(F2F_h, calculatedF);
          if(fgreater(calculatedF, bound)){
            UpdateNextBound(bound, calculatedF);
            ignoreList.push_back(childID);
            heuristicList.erase(heuristicList.begin()+i);
            i--;
          }
        }
      }
      if(heuristicList.size() ==0){
        for (uint64_t childID : ignoreList){
          heuristicList.push_back(childID);
        }
        ignoreList.clear();
        return F2F_h;
      }
    }
    else{
      for (AStarOpenClosedDataWithF<state> astarState : statesList.getElements()){
        if(astarState.g == perimeterG){ //perimeter frontier should be fixed and this should be updated to <= instead of ==
          F2F_h = std::min(F2F_h,heuristic->HCost(currState, astarState.data) + astarState.g);
        }
      }
      if (fgreater(g+F2F_h, bound)){
        UpdateNextBound(bound, g+F2F_h);
        //printf("Stopping at (%d, %d). g=%f h=%f\n", currState>>16, currState&0xFFFF, g, h);
        return F2F_h;
      }
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
		if (neighbors[x] == parent || (isDuplicateDetection && statesList.Lookup(env->GetStateHash(neighbors[x]), childID) == kClosedList)) {
			continue;
		}
		thePath.push_back(neighbors[x]);
		double edgeCost = env->GCost(currState, neighbors[x]);
		double childH = DoIteration(env, currState, neighbors[x], thePath, bound,
																g+edgeCost, maxH - edgeCost);
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
	for (uint64_t childID : ignoreList){
		heuristicList.push_back(childID);
	}
	ignoreList.clear();
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


template <class state, class action, bool verbose>
bool IDAStar<state, action, verbose>::SFBDS(SearchEnvironment<state, action> *env, state from, state to,
							 std::vector<state> &thePath, int secondsLimit, int lookahead, bool isAlternating, bool storeLookup)
{
	if(verbose){
		printf("\t\tStarting to solve with SFBDS\n");
	}
	if (!storedHeuristic)
		heuristic = env;
	nextBound = 0;
	nodesExpanded = nodesTouched = dAstarLastIterExpansions = 0;
	thePath.resize(0);
  double initialHeuristic = heuristic->HCost(from, to);
	UpdateNextBound(0, initialHeuristic);
	goal = to;
	start = from;
  std::vector<state> startPath;
  std::vector<state> goalPath;
	startPath.push_back(from);
  goalPath.push_back(to);
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
		if (verbose)
			printf("Nodes expanded: %llu(%llu)\n", nodesExpanded-nodesExpandedSoFar, nodesExpanded);
		if (DoSFBDSIteration(env, from,from, to,to, startPath,goalPath, nextBound, 0, 0,lookahead,1,initialHeuristic,isAlternating,storeLookup,true)){
			if (verbose){
				printf("Nodes expanded: %llu(%llu)\n", nodesExpanded-nodesExpandedSoFar, nodesExpanded);
			}
			necessaryExpansions += nodesExpandedSoFar;
			break;
		}
    if (verbose){
			printf("Nodes expanded: %llu(%llu)\n", nodesExpanded-nodesExpandedSoFar, nodesExpanded);
		}
		previousIterationExpansions = nodesExpanded-nodesExpandedSoFar;
		nodesExpandedSoFar = nodesExpanded;
		PrintGHistogram();
	}
	dAstarExpansions = previousIterationExpansions + dAstarLastIterExpansions;
  goalPath.pop_back();
  std::reverse(goalPath.begin(),goalPath.end());
	startPath.insert( startPath.end(), goalPath.begin(), goalPath.end() );
  thePath = startPath;
	return true;
}

template <class state, class action, bool verbose>
bool IDAStar<state, action, verbose>::DoSFBDSIteration(SearchEnvironment<state, action> *env,
										   state parentStart, state currStart, state parentGoal, state currGoal,
										   std::vector<state> &startPath, std::vector<state> &goalPath, double bound, double g_from, double g_to, int k,uint64_t b,double h, bool isAlternating,bool storeLookup, bool isForward)
{
	if (currStart == currGoal && !fgreater(g_from + g_to,bound)){
		return true;
  }
  //double h = heuristic->HCost(currStart, currGoal);
	// path max

	if (fgreater(g_from +h +g_to, bound))
	{
		UpdateNextBound(bound, g_from +h +g_to);
		return false;
	}
  std::vector<state> forwardFrontierStates;
  std::vector<double> forwardFrontierH;
  std::vector<state> backwardFrontierStates;
  std::vector<double> backwardFrontierH;
  bool isSolved = false;
 if (!isAlternating){
    uint64_t forwardValue = computeKlook(env,parentStart,currStart,currGoal,bound,g_from,g_to,g_from +h +g_to,k,startPath,isSolved,b,storeLookup,forwardFrontierStates,forwardFrontierH);
    if (isSolved){
      return true;
    }
    uint64_t backwardValue = computeKlook(env,parentGoal, currGoal, currStart,bound,g_to,g_from,g_from +h +g_to,k,goalPath,isSolved,b,storeLookup,backwardFrontierStates,backwardFrontierH);
    if (isSolved){
      return true;
    }
    isForward = forwardValue <= backwardValue;
 }
  if (isForward){
    std::vector<state> startNeighbors;
    if (storeLookup){
      startNeighbors = forwardFrontierStates;
    }
    else{
      env->GetSuccessors(currStart, startNeighbors);
      nodesTouched += startNeighbors.size();
      nodesExpanded++;
    }

    for (unsigned int x = 0; x < startNeighbors.size(); x++){
      uint64_t childID;
      if (startNeighbors[x] == parentStart) {
        continue;
      }
      double nextH;
      if (storeLookup){
        nextH = forwardFrontierH[x];
      }
      else{
        nextH = heuristic->HCost(startNeighbors[x], currGoal);
      }
      startPath.push_back(startNeighbors[x]);
      double edgeCost = env->GCost(currStart, startNeighbors[x]);
      isSolved = DoSFBDSIteration(env, currStart, startNeighbors[x], parentGoal,currGoal,startPath,goalPath, bound,
                                  g_from+edgeCost, g_to,k,startNeighbors.size(),nextH,isAlternating,storeLookup,!isForward);
      if(isSolved){
        return true;
      }
      startPath.pop_back();
    }
    return false;
  }
  else{
    std::vector<state> goalNeighbors;
    if (storeLookup){
      goalNeighbors = backwardFrontierStates;
    }
    else{
      env->GetSuccessors(currGoal, goalNeighbors);
      nodesTouched += goalNeighbors.size();
      nodesExpanded++;
    }
    for (unsigned int x = 0; x < goalNeighbors.size(); x++){
      uint64_t childID;
      if (goalNeighbors[x] == parentGoal) {
        continue;
      }
      double nextH;
      if (storeLookup){
        nextH = backwardFrontierH[x];
      }
      else{
        nextH = heuristic->HCost(currStart, goalNeighbors[x]);
      }
      goalPath.push_back(goalNeighbors[x]);
      double edgeCost = env->GCost(currGoal, goalNeighbors[x]);
      isSolved = DoSFBDSIteration(env, parentStart,currStart, currGoal,goalNeighbors[x],startPath,goalPath, bound,
                                  g_from,g_to+edgeCost,k,goalNeighbors.size(),nextH,isAlternating,storeLookup,!isForward);
      if(isSolved){
        return true;
      }
      goalPath.pop_back();
    }
    return false;
  }
  return false;
}

template <class state, class action, bool verbose>
double IDAStar<state, action, verbose>::computeKlook(SearchEnvironment<state, action> *env,state parent, state curr, state opposite,double bound,double g_curr,double g_opposite,double prevf,int k,std::vector<state> &thePath,bool &isSolved, uint64_t b, bool storeLookup, std::vector<state> &stateStorage, std::vector<double> &hStorage){
  
  if (curr == opposite && !fgreater(g_curr + g_opposite,bound)){
		isSolved = true;
    return 0;
  }
  double h = heuristic->HCost(curr, opposite);
	// path max
  double nodesBound = g_curr + g_opposite +h;
	if (fgreater(nodesBound, bound))
	{
    UpdateNextBound(bound, g_curr + g_opposite +h);
		return 0;
	}
  uint64_t value = 0;
  if (nodesBound == bound){
    value = 1;
  }
  else if (nodesBound == prevf){
    value = b;
  }
  else{
    value = 1;
  }
  if (k==0){
    if (storeLookup){
      stateStorage.push_back(curr);
      hStorage.push_back(h);
    }
    return value;
  }
  std::vector<state> neighbors;
  env->GetSuccessors(curr, neighbors);
  nodesTouched += neighbors.size();
  nodesExpanded++;
  for (unsigned int x = 0; x < neighbors.size(); x++){
    uint64_t childID;
    if (neighbors[x] == parent) {
      continue;
    }

    thePath.push_back(neighbors[x]);
    double edgeCost = env->GCost(curr, neighbors[x]);
    value += computeKlook(env, curr,neighbors[x],opposite, bound,
                                g_curr+edgeCost,g_opposite,nodesBound,k-1,thePath,isSolved,b, storeLookup, stateStorage, hStorage);
    if(isSolved){
      return value;
    }
    thePath.pop_back();
  }
  return value;
}


#endif

