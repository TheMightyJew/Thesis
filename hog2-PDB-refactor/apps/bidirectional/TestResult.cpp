#include "TestResult.h"

TestResult::TestResult(TestInfo testInfo) : m_testInfo(testInfo), m_initialHeuristic(0), m_algorithmInfo(""), m_solutionCost(-1), m_statesExpanded(0), m_neccessaryStatesExpanded(0), m_maxStatesInMemory(0), m_iterations(0), m_timeElapsed(0) {}

string TestResult::csvSerializeHeaders()
{
    return TestInfo::csvSerializeHeaders() + ", Initial Heuristic, Algorithm Info, Solution Length, States Expanded, Neccessary States Expanded, Max States In Memory, Iterations, Time Elapsed";
}

string TestResult::csvSerialize()
{
    return m_testInfo.csvSerialize() + ",\"" + to_string(m_initialHeuristic) + "\",\"" + m_algorithmInfo + "\",\"" + to_string(m_solutionCost) + "\",\"" + to_string(m_statesExpanded) + "\",\"" + to_string(m_neccessaryStatesExpanded) + "\",\"" + to_string(m_maxStatesInMemory) + "\",\"" + to_string(m_iterations) + "\",\"" + to_string(m_timeElapsed) + "\"";
}
