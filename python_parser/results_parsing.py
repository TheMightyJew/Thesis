import pandas as pd
import numpy as np
import sys
import os

analysis_dir = 'Analysis'
res_dir = 'PancakeSorting'
path = analysis_dir + '/' + res_dir + '/'
filename = ''
try:
    arg = sys.argv[1]
    path = os.path.dirname(arg) + '/'
    filename = os.path.basename(arg)
    filenames = [filename]
except:
    filenames = ['three_maps', 'one_map']
for filename in filenames:
    filename = str(filename)
    fileName = filename.replace('.txt', '')
    file = open(path + fileName + ".txt", "r")
    if res_dir == 'PancakeSorting':
        cols = ['Number Of Pancakes', 'Gap']
        errorCols = ['Number Of Pancakes', 'Gap']
        analysisCols = ['Gap']
    else:
        cols = []
        errorCols = []
        analysisCols = []
    if res_dir != 'Grid':
        cols += ['Start state', 'Goal state', 'Initial Heuristic']
    else:
        cols += ['Map Name']
    cols += ['Problem ID', 'Algorithm',
             'Memory', 'Status', 'States Expanded', 'Necessary Expansions', 'Iterations', 'Runtime(seconds)']
    errorCols += ['Problem ID', 'Algorithm', 'Memory', 'Status', 'solLength',
                  'realSolLength']
    analysisCols += ['Algorithm', 'Failed(Out Of 100)', 'Mean Expansions', 'Mean Necessary Expansions',
                     'Mean Iterations', 'Mean Runtime(Seconds)']
    resultsDF = pd.DataFrame()
    errorsDF = pd.DataFrame()
    resDict = dict.fromkeys(cols)
    errorsDict = dict.fromkeys(errorCols)
    gapStr = 'GAP-0'
    map_name = ''
    resDict['Problem ID'] = 0
    errorsDict['Problem ID'] = 0
    necessary = 0
    for line in file:
        isError = False
        line = line.replace('\t', '').replace('\n', '')
        splittedLine = line.replace(gapStr, ' ' + gapStr + ' ').split()
        skipList = ['command line arguments', 'U2 F2 U D R2 U- D- F2', './bidirectional']
        if 'command line arguments' in line or 'U2 F2 U D R2 U- D- F2' in line or './bidirectional' in line:
            continue
        if 'Test' in line:
            if 'TestPancake' in line:
                resDict['Number Of Pancakes'] = int(
                    splittedLine[splittedLine.index('TestPancake:(Pancakes:') + 1].replace(',', ''))
                resDict['Gap'] = int(splittedLine[splittedLine.index('Gap:') + 1].replace(')', '').replace(',', ''))
                gapStr = 'GAP-' + str(resDict['Gap'])
        elif 'Loading' in line and res_dir == 'Grid':
            a = line[:line.index('.map')].split('/')[-1]
            map_name = line
        elif 'Problem' in line:
            problemID = int(splittedLine[splittedLine.index('Problem') + 1])
            MMsolLength = 0
            AsolLength = 0
            solLength = 0
            realSolLength = 0
            resDict['Problem ID'] = problemID
            if res_dir != 'Grid':
                line = next(file).replace('\t', '').replace('\n', '')
                resDict['Start state'] = line[line.index(':') + 1:-1]
                line = next(file).replace('\t', '').replace('\n', '')
                resDict['Goal state'] = line[line.index(':') + 1:-1]
                line = next(file).replace('\t', '').replace('\n', '')
                resDict['Initial Heuristic'] = float(line[line.index('heuristic') + len('heuristic') + 1:])
        else:
            if line[0] == line[-1] == '_' or 'completed' in line:
                continue
            resDict['Algorithm'] = splittedLine[0]
            if 'necessary' in line:
                necessary = int(splittedLine[splittedLine.index('necessary;') - 1])
            resDict['Necessary Expansions'] = necessary
            if 'MBBDS' in resDict['Algorithm']:
                if 'iterations;' in splittedLine:
                    resDict['Iterations'] = int(splittedLine[splittedLine.index('iterations;') - 1])
                elif 'iterations' in splittedLine:
                    resDict['Iterations'] = int(splittedLine[splittedLine.index('iterations') - 1])
            else:
                resDict['Iterations'] = 0
            if line.startswith('MM found'):
                MMsolLength = float(splittedLine[splittedLine.index('length') + 1].replace(';', ''))
            if line.startswith('A* found'):
                AsolLength = float(splittedLine[splittedLine.index('length') + 1].replace(';', ''))
            if 'memory' in line:
                if 'MBBDS' in resDict['Algorithm'] or '+' in resDict['Algorithm'] or 'BAI' in resDict[
                    'Algorithm'] or 'IDTHSpTrans' in resDict['Algorithm'] or 'BFBDS' in resDict['Algorithm']:
                    memoryStr = 'Memory_Percentage='
                    resDict['Memory'] = int(splittedLine[splittedLine.index('memory') + 2])
                    resDict['Algorithm'] += '(' + line[line.index(memoryStr) + len(memoryStr):][:4] + ')'
                else:
                    try:
                        resDict['Memory'] = float(splittedLine[splittedLine.index('using') + 1])
                    except:
                        a = 1
            algos2Check = ['MBBDS', 'A*+IDA*(', 'A*+IDA*_Reverse(', 'MM+IDMM', 'A*+IDMM', 'IDMM', 'BAI', 'IDTHSpTrans ',
                           'IDTHSpTrans_NDD']
            # if 'length' in line and not line.startswith('MM ') and not line.startswith('A* '):
            if 'length' in line and not line.startswith('A* '):
                solLength = float(splittedLine[splittedLine.index('length') + 1].replace(';', ''))
                # realSolLength = max(AsolLength, MMsolLength)
                realSolLength = min(AsolLength, MMsolLength)
                # if solLength not in [AsolLength, MMsolLength]:
                if solLength != AsolLength:
                    isError = True
            else:
                resDict['Memory'] = '-'
            if 'fail' in line:
                resDict['Status'] = 0
                resDict['States Expanded'] = '-'
                resDict['Runtime(seconds)'] = '-'
            else:
                resDict['Status'] = 1
                try:
                    resDict['States Expanded'] = int(splittedLine[splittedLine.index('expanded;') - 1])
                except:
                    print('error:', line)
                if 'elapsed' in splittedLine:
                    resDict['Runtime(seconds)'] = float(
                        splittedLine[splittedLine.index('elapsed') - 1].replace('s', ''))
                elif 'elapsed;' in splittedLine:
                    resDict['Runtime(seconds)'] = float(
                        splittedLine[splittedLine.index('elapsed;') - 1].replace('s', ''))
            resultsDF = resultsDF.append(resDict, ignore_index=True)
            if isError:
                algoName = resDict['Algorithm']
                if '(' in algoName:
                    algoName = algoName[:algoName.index('(')]
                if res_dir == 'PancakeSorting':
                    errorsDict = {'Number Of Pancakes': resDict['Number Of Pancakes'], 'Gap': resDict['Gap'],
                                  'Problem ID': resDict['Problem ID'], 'Algorithm': algoName,
                                  'Memory': resDict['Memory'],
                                  'Status': resDict['Status'],
                                  'solLength': solLength,
                                  'realSolLength': realSolLength}
                elif res_dir == 'STP':
                    errorsDict = {'Problem ID': resDict['Problem ID'], 'Algorithm': algoName,
                                  'Memory': resDict['Memory'],
                                  'Status': resDict['Status'],
                                  'solLength': solLength,
                                  'realSolLength': realSolLength}
                elif res_dir == 'Grid':
                    errorsDict = {'Map Name': map_name, 'Problem ID': resDict['Problem ID'], 'Algorithm': algoName,
                                  'Memory': resDict['Memory'],
                                  'Status': resDict['Status'],
                                  'solLength': solLength,
                                  'realSolLength': realSolLength}
                errorsDF = errorsDF.append(errorsDict, ignore_index=True)

    resultsDF = resultsDF[cols]
    resultsDF.to_csv(path + fileName + '_results.csv')
    errorsDF.drop_duplicates(inplace=True)
    if len(errorsDF) > 0:
        if res_dir == 'PancakeSorting':
            errorsDF.sort_values(['Number Of Pancakes', 'Gap', 'Problem ID', 'Algorithm', 'Memory'], ascending=True,
                                 inplace=True)
        else:
            errorsDF.sort_values(['Problem ID', 'Algorithm', 'Memory'], ascending=True,
                                 inplace=True)
        errorsDF = errorsDF[errorCols]
    errorsDF.to_csv(path + fileName + '_errors.csv')
    analysisDF = pd.DataFrame()
    analysisDict = dict.fromkeys(analysisCols)
    MBBDSFailSum = None
    failedRows = (resultsDF['Status'] == 0)
    if res_dir == 'PancakeSorting':
        for gap in resultsDF['Gap'].unique():
            analysisDict['Gap'] = gap
            gapRows = (resultsDF['Gap'] == gap)
            finished_problems_Rows = (resultsDF['Problem ID'].isin(
                set(resultsDF[gapRows]['Problem ID']) - set(resultsDF[gapRows & (failedRows)]['Problem ID'])))
            for algo in resultsDF['Algorithm'].unique():
                analysisDict['Algorithm'] = algo
                algoRows = (resultsDF['Algorithm'] == algo)
                if MBBDSFailSum == algo:
                    analysisDict['Failed(Out Of 100)'] += len(resultsDF[gapRows & algoRows & failedRows])
                else:
                    analysisDict['Failed(Out Of 100)'] = len(resultsDF[gapRows & algoRows & failedRows])
                analysisDict['Mean Expansions'] = int(
                    np.mean(resultsDF[gapRows & algoRows & finished_problems_Rows]['States Expanded']))
                analysisDict['Mean Necessary Expansions'] = int(
                    np.mean(resultsDF[gapRows & algoRows & finished_problems_Rows]['Necessary Expansions']))
                analysisDict['Mean Iterations'] = int(
                    np.mean(resultsDF[gapRows & algoRows & finished_problems_Rows]['Iterations']))
                analysisDict['Mean Runtime(Seconds)'] = round(
                    np.mean(resultsDF[gapRows & algoRows & finished_problems_Rows]['Runtime(seconds)']), 5)
                analysisDF = analysisDF.append(analysisDict, ignore_index=True)
                MBBDSFailSum = algo
    else:
        finished_problems_Rows = (resultsDF['Problem ID'].isin(
            set(resultsDF['Problem ID']) - set(resultsDF[failedRows]['Problem ID'])))
        for algo in resultsDF['Algorithm'].unique():
            analysisDict['Algorithm'] = algo
            algoRows = (resultsDF['Algorithm'] == algo)
            if MBBDSFailSum == algo:
                analysisDict['Failed(Out Of 100)'] += len(resultsDF[algoRows & failedRows])
            else:
                analysisDict['Failed(Out Of 100)'] = len(resultsDF[algoRows & failedRows])
            analysisDict['Mean Expansions'] = int(
                np.mean(resultsDF[algoRows & finished_problems_Rows]['States Expanded']))
            analysisDict['Mean Necessary Expansions'] = int(
                np.mean(resultsDF[algoRows & finished_problems_Rows]['Necessary Expansions']))
            analysisDict['Mean Iterations'] = int(
                np.mean(resultsDF[algoRows & finished_problems_Rows]['Iterations']))
            analysisDict['Mean Runtime(Seconds)'] = round(
                np.mean(resultsDF[algoRows & finished_problems_Rows]['Runtime(seconds)']), 5)
            analysisDF = analysisDF.append(analysisDict, ignore_index=True)
            MBBDSFailSum = algo

    analysisDF[analysisCols].to_csv(path + fileName + '_analysis.csv')