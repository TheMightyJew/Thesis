/*
 *  IDTHSwTrans.h
 *
 *  Created by Steven Danishevski on MARCH 20.
 *
 */

#ifndef IDTHS_TRANS_H
#define IDTHS_TRANS_H

#include <iostream>
#include <unordered_set>
#include "SearchEnvironment.h"
#include <math.h>


template<typename state>
class transpostionNode {
public:
	transpostionNode(const state &theData, double gCost, double hCost, double error, uint64_t hashValue, bool open):data(theData), error(error), g(gCost), h(hCost), hash(hashValue), isOpen(open) {}
	state data;
	double g;
	double h;
  double error;
  uint64_t hash;
  bool isOpen;
};

template<typename state, class dataStructure = transpostionNode<state> >
class TranspositionTable {
public:
  TranspositionTable() {}
	TranspositionTable(bool hash) : useHash(hash){}
	~TranspositionTable() {}
	uint64_t AddNode(const state &val, uint64_t hash, double g, double h, double error);
	bool Lookup(uint64_t hashKey, uint64_t &objKey) const;
	inline dataStructure &Lookup(uint64_t objKey) { return elements[objKey]; }
  inline dataStructure &LookupOpen(uint64_t objKey) { return elements[currentOpenSet[objKey]]; }
  inline uint64_t &GetOpenKey(uint64_t objKey) { return currentOpenSet[objKey]; }
	inline const dataStructure &Lookat(uint64_t objKey) const { return elements[objKey]; }
  inline void Reset(){ elements.clear(); table.clear(); currentOpenSet.clear();}
  inline void close(){ elements.clear(); table.clear(); currentOpenSet.clear();}
  inline void reOpenNode(uint64_t key){ elements[key].isOpen = true; currentOpenSet.push_back(key);}
  typedef __gnu_cxx::hash_map<uint64_t, uint64_t, AHash64> IndexTable;
	IndexTable table;
  uint64_t size() {return elements.size();}
  uint64_t GetNumOpenElements() {return currentOpenSet.size();}
  void closeNodeAtIndex(uint64_t index);
  void closeNodeAtElementIndex(uint64_t elementIndex);
  
private:
  std::vector<dataStructure> elements;
  std::vector<uint64_t> currentOpenSet;
  bool useHash;
};


template<typename state, class dataStructure>
uint64_t TranspositionTable<state, dataStructure>::AddNode(const state &val, uint64_t hash, double g, double h, double error)
{
	elements.push_back(dataStructure(val, g, h, error,hash,true));
  if (useHash){
    table[hash] = elements.size()-1; // hashing to element list location
  }
  currentOpenSet.push_back(elements.size()-1);
	return elements.size()-1;
}
template<typename state, class dataStructure>
void TranspositionTable<state, dataStructure>::closeNodeAtIndex(uint64_t index)
{
  elements[currentOpenSet[index]].isOpen = false;
	currentOpenSet.erase(currentOpenSet.begin()+index);
}

template<typename state, class dataStructure>
void TranspositionTable<state, dataStructure>::closeNodeAtElementIndex(uint64_t elementIndex)
{
  elements[elementIndex].isOpen = false;
  currentOpenSet.erase(std::remove(currentOpenSet.begin(),currentOpenSet.end(), elementIndex), currentOpenSet.end());
}

template<typename state, class dataStructure>
bool TranspositionTable<state, dataStructure>::Lookup(uint64_t hashKey, uint64_t &objKey) const
{
	typename IndexTable::const_iterator it;
	it = table.find(hashKey);
	if (it != table.end())
	{
		objKey = (*it).second;
		return true;
	}
	return false;
}

template <class state, class action, bool verbose = true, class table = TranspositionTable<state>>
class IDTHSwTrans {
public:
	IDTHSwTrans(bool front2frontH=false, bool isConsistent = false, bool isUpdateByWorkload = false, double smallestEdge=1, bool useHash = true) { this->front2frontH = front2frontH; this->isConsistent = isConsistent; this->smallestEdge = smallestEdge; this->isUpdateByWorkload = isUpdateByWorkload; this->useHash = useHash;}
	virtual ~IDTHSwTrans() {}
	bool GetPath(SearchEnvironment<state, action>* env, state fromState, state toState, /*std::vector<state> &thePath,*/ int secondsLimit=600, unsigned long availableStorage = 0);
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
  bool useHash;
	bool DoIterationForward(SearchEnvironment<state, action>* env,
	state parent, state currState, uint64_t hashValue, double g, uint64_t &midState);//, std::vector<uint64_t> &currPath, std::vector<state> &backwardPath);	 
	bool DoIterationBackward(SearchEnvironment<state, action>* env, state parent, state currState, double g, uint64_t& midState);//, std::vector<state>& backwardPath);
	double updateBoundByFraction(double boundToSplit,double p = 0.5, bool isInteger = true);
	double updateBoundByWorkload(double newbound, double prevBound, double oldForwardBound,uint64_t forwardLoad,uint64_t backwardLoad);
	double nextBound[2];
	int const forwardLoc = 0;
	int const backwardLoc = 1;
	void UpdateNextBound(double currBound, double fCost, int loc);
	state originGoal;
	state originStart;
	unsigned long dMMLastIterExpansions = 0;
	unsigned long necessaryExpansions = 0;
	unsigned long nodesExpandedSoFar = 0;
	unsigned long previousIterationExpansions = 0;
	uint64_t forwardExpandedInLastIter = 0;
	uint64_t backwardExpandedInLastIter = 0;
  
  unsigned long availableStorage = 0;
	

  table transTable;
	bool readyOpenLists = false;
	bool front2frontH;
	bool firstBounds;
	bool isConsistent;
	double smallestEdge;
	bool isUpdateByWorkload;
	std::chrono::steady_clock::time_point startTimeTest;

};

template <class state, class action, bool verbose, class table>
double IDTHSwTrans<state, action, verbose,table>::updateBoundByWorkload(double newbound, double prevBound, double oldForwardBound, uint64_t forwardLoad,uint64_t backwardLoad){
	//printf("forwardLoad: %llu,backwardLoad: %llu,oldForwardBound%d\n",forwardLoad,backwardLoad,oldForwardBound);
  //printf("update3: %1.0f|%1.0f|%1.0f|%llu|%llu\n", newbound, prevBound, oldForwardBound, forwardLoad, backwardLoad);
	if (forwardLoad <= backwardLoad){
		return oldForwardBound + newbound - prevBound;
	}
	else{
		return oldForwardBound;
	}
}

template <class state, class action, bool verbose, class table>
bool IDTHSwTrans<state, action, verbose, table>::GetPath(SearchEnvironment<state, action>* env,
state fromState, state toState, /*std::vector<state> &thePath,*/ int secondsLimit, unsigned long availableStorage)
{
  //printf("start\n");
  //thePath.clear();
  this->availableStorage = availableStorage;
	auto startTime = std::chrono::steady_clock::now();
	nodesExpanded = nodesTouched = 0;
	originStart = fromState;
	originGoal = toState;
	double initialHeuristic = env->HCost(fromState, toState);
	fBound = 	nextBound[forwardLoc] = nextBound[backwardLoc] = std::max(smallestEdge,initialHeuristic);
	forwardBound = ceil(fBound / 2) - smallestEdge;
	unsigned long nodesExpandedSoFar = 0;
	unsigned long previousIterationExpansions = 0;
  uint64_t midState;
  //std::vector<uint64_t> currPath;
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
    uint64_t hash;
    if (useHash){
      hash = env->GetStateHash(originStart);
    }
    //printf("bound=%1.0f, forwardBound%1.0f\n", fBound, forwardBound);
		bool solved = DoIterationForward(env, originStart, originStart, hash, 0, midState);//, currPath, thePath);
		if(verbose){
			printf("Nodes expanded: %d(%d)\n", nodesExpanded-nodesExpandedSoFar, nodesExpanded);
		}

    if (!solved && transTable.GetNumOpenElements() > 0){
      solved = DoIterationBackward(env, originGoal, originGoal, 0, midState);//, thePath);
    }
    if (solved) {
			dMMExpansions = previousIterationExpansions + dMMLastIterExpansions;
			necessaryExpansions += nodesExpandedSoFar;
      /*
      uint64_t oldMidState =  std::numeric_limits<uint64_t>::max();
      while (midState != oldMidState){
        oldMidState = midState;
        thePath.push_back(transTable.Lookup(midState).data);
        midState = transTable.Lookup(midState).parentID;
      }
			std::reverse(thePath.begin(),thePath.end());
      */
      return true;
		}
		else{
      //printf("%1.0f|%1.0f|%1.0f|%1.0f|%llu|%llu\n", forwardBound, fBound, nextBound[forwardLoc], nextBound[backwardLoc], forwardExpandedInLastIter, backwardExpandedInLastIter);
      transTable.Reset();
			double nextFbound = std::max(nextBound[forwardLoc],nextBound[backwardLoc]);
      nextBound[forwardLoc] = nextBound[backwardLoc] = nextFbound;
			if (!isUpdateByWorkload){
				forwardBound = ceil(nextFbound / 2) - smallestEdge;
			} 
			else{
				forwardBound = updateBoundByWorkload(nextFbound, fBound, forwardBound, forwardExpandedInLastIter, backwardExpandedInLastIter);
			}
			fBound = nextFbound;
			forwardExpandedInLastIter = 0;
			backwardExpandedInLastIter = 0;
		}
		previousIterationExpansions = nodesExpanded-nodesExpandedSoFar;
		nodesExpandedSoFar = nodesExpanded;
	}
	return false;
}


template <class state, class action, bool verbose, class table>
bool IDTHSwTrans<state, action, verbose, table>::DoIterationForward(SearchEnvironment<state, action>* env,
	state parent, state currState, uint64_t hashValue, double g, uint64_t &midState)//, std::vector<uint64_t> &currPath, std::vector<state> &backwardPath)
{
	double h = env->HCost(currState, originGoal);
	if (fgreater(g + h, fBound)){
		UpdateNextBound(fBound, g + h,forwardLoc);
		return false;
	}

	if (g > forwardBound) {
		//backwardBound = fBound - g - smallestEdge;
    double error = 0;
		if (isConsistent){
		  error = g - env->HCost(currState,originStart);
		}
    uint64_t frontierKey = transTable.AddNode(currState,hashValue, g, h, error);
    if (transTable.size() > availableStorage){
      bool isSolutionFound = DoIterationBackward(env, originGoal, originGoal, 0, midState);//,backwardPath);
      if (isSolutionFound){
        return true;
      }
      else{
        transTable.Reset();
        return false;
      }
    }
    return false;	
	}
	std::vector<state> neighbors;
	env->GetSuccessors(currState, neighbors);
	nodesTouched += neighbors.size();
	nodesExpanded++;
  forwardExpandedInLastIter++;
	if(g + h == fBound){
		dMMLastIterExpansions++;
	}
	for (unsigned int x = 0; x < neighbors.size(); x++){
		uint64_t childID;
		double edgeCost = env->GCost(currState, neighbors[x]);
		if (neighbors[x] == parent){
      continue;
    }
    uint64_t hash;
    if (useHash){
      hash = env->GetStateHash(neighbors[x]);
    }
    if (useHash && transTable.Lookup(hash, childID)){
      if (!transTable.Lookup(childID).isOpen){
        continue;
      }
      if (transTable.Lookup(childID).g <= g + edgeCost) {
        continue;
      }
      else{
        transTable.Lookup(childID).g = g + edgeCost;
        if (g + edgeCost <= forwardBound)
        transTable.closeNodeAtElementIndex(childID);
      }
		}
		if (DoIterationForward(env, currState, neighbors[x],hash, g + edgeCost, midState)){//, currPath, backwardPath)) {
			return true;
		}
	}
  //currPath.pop_back();
	return false;
}

template <class state, class action, bool verbose, class table>
bool IDTHSwTrans<state, action, verbose, table>::DoIterationBackward(SearchEnvironment<state, action>* env,
	state parent, state currState, double g, uint64_t& midState)//, std::vector<state>& backwardPath)
{	

  double fPossibleBound = 0;
  //backwardPath.push_back(currState);
  /*
    if (front2frontH){
    double h = env->HCost(currState, possibleMidState);
    fPossibleBound = std::max(fPossibleBound, g + h + possibleMidStateG);
    if (fgreater(fPossibleBound, fBound) || g > backwardBound){
      UpdateNextBound(fBound, fPossibleBound,backwardLoc);
      return false;
    }
  }
  */

  std::vector<uint64_t> ignoreList;
  for (uint64_t i = transTable.GetNumOpenElements()-1; i != static_cast<uint64_t>(-1);i--){
    state &possibleMidState = transTable.LookupOpen(i).data;
    double possibleMidStateG = transTable.LookupOpen(i).g;
   // double otherError = transTable.LookupOpen(i).error;
    //double otherH = transTable.LookupOpen(i).h;
    backwardBound = fBound - possibleMidStateG - smallestEdge;
    if (currState == possibleMidState && !fgreater(g + possibleMidStateG,fBound)) {
      pathLength = g + possibleMidStateG;
      midState = i;
      return true;
    }
    double computedF = g + possibleMidStateG + env->HCost(currState,possibleMidState);
    if(fgreater(computedF,fBound) || g > backwardBound){
      ignoreList.push_back(transTable.GetOpenKey(i));
      transTable.closeNodeAtIndex(i);
      UpdateNextBound(fBound, computedF,backwardLoc);
      continue;
    }
    /*
    double originalH = env->HCost(currState, originStart);
    //fPossibleBound = std::max(g + originalH + otherError, fPossibleBound);
    fPossibleBound = std::max(g + originalH, fPossibleBound);
    if (fgreater(fPossibleBound, fBound) || g > backwardBound){
      UpdateNextBound(fBound, fPossibleBound, backwardLoc);
      ignoreList.push_back(transTable.GetOpenKey(i));
      transTable.closeNodeAtIndex(i);
      i--;
      continue;
    }
    
    if (isConsistent){
      double error = g - env->HCost(currState,originGoal);
      //fPossibleBound = std::max(possibleMidStateG + otherH + error, fPossibleBound); 
      fPossibleBound = std::max(possibleMidStateG + otherH, fPossibleBound); 
      if (fgreater(fPossibleBound, fBound) || g > backwardBound){
        UpdateNextBound(fBound, fPossibleBound,backwardLoc);
        ignoreList.push_back(transTable.GetOpenKey(i));
        transTable.closeNodeAtIndex(i);
        i--;
        continue;
      }
    }
*/    
  }
  if(transTable.GetNumOpenElements() ==0){
    for (uint64_t i : ignoreList){
      transTable.reOpenNode(i);
    }
    ignoreList.clear();
    //backwardPath.pop_back();
    return false;
  }
   
  std::vector<state> neighbors;
  env->GetSuccessors(currState, neighbors);
  nodesTouched += neighbors.size();
  nodesExpanded++;
  backwardExpandedInLastIter++;
  if(fPossibleBound == fBound){
    dMMLastIterExpansions++;
  }
  for (unsigned int x = 0; x < neighbors.size(); x++){
    uint64_t childID;
    double edgeCost = env->GCost(currState, neighbors[x]);
    if (neighbors[x] == parent) {
      continue;
    }
    if (DoIterationBackward(env, currState, neighbors[x], g + edgeCost, midState)){//, backwardPath)) {
      return true;
    }
  }
  for (uint64_t i : ignoreList){
    transTable.reOpenNode(i);
  }
  ignoreList.clear();
  //backwardPath.pop_back();
  return false;
}

template <class state, class action, bool verbose, class table>
void IDTHSwTrans<state, action, verbose, table>::UpdateNextBound(double currBound, double fCost, int loc)
{
	if (!fgreater(nextBound[loc], currBound))
	{
		nextBound[loc] = std::max(fCost,currBound);
	}
	else if (fgreater(fCost, currBound) && fless(fCost, nextBound[loc]))
	{
		nextBound[loc] = fCost;
	}
}



#endif