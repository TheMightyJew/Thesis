#ifndef TestResult_H
#define TestResult_H
 
#include <iostream>
#include <string>
#include "TestInfo.h"
#include "CsvSerializable.h"

using namespace std;

class TestResult
{
public:
    TestInfo m_testInfo;
    double m_initialHeuristic;
    string m_algorithmInfo;
    double m_solutionCost;
    unsigned int m_statesExpanded;
    unsigned int m_neccessaryStatesExpanded;
    unsigned int m_maxStatesInMemory;
    unsigned int m_iterations;
    float m_timeElapsed;

    TestResult(TestInfo testInfo);
    static string csvSerializeHeaders();
    string csvSerialize();
};
 
#endif