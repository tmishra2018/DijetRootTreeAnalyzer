import os

command = 'bash command_HTcondor.sh'

configs={}
cutfiles={}

samples = {
    'B' : 'Run2017B-09Aug2019_UL2017-v1',
    'C' : 'Run2017C-09Aug2019_UL2017-v1',
    'D' : 'Run2017D-09Aug2019_UL2017-v1',
    'E' : 'Run2017E-09Aug2019_UL2017-v1',
    'F' : 'Run2017F-09Aug2019_UL2017-v1',
    'MC1': 'GJets_HT-40To100_RunIISummer19MiniAOD-106X',
    'MC2': 'GJets_HT-100To200_RunIISummer19MiniAOD-106X',
    'MC3': 'GJets_HT-200To400_RunIISummer19MiniAOD-106X',
    'MC4': 'GJets_HT-400To600_RunIISummer19MiniAOD-106X',
    'MC5': 'GJets_HT-600ToInf_RunIISummer19MiniAOD-106X',
}

tags = {
    'withoutL2Res' : 'wo_L2Res', 
    'onlyL2Res' : 'only_L2Res', 
    'L2L3Res' : 'L2L3Res', 
    'JER' : 'JER'
}

JERCs = tags.keys()
runs = samples.keys()

for run in runs:
    for JERC in JERCs:
        key = '_'.join([run,JERC])
        if 'MC' in run:
            cutfiles[key] = cutfiles['_'.join(['B',JERC])]
        else:
            cutfiles[key] = 'config/cutFile_Run2017UL_{run}_{JERC}.txt'.format(
            run=run,
            JERC=JERC
        )
        configs[key] = '{cutfile} {dirlist} {tag}'.format(
            cutfile = cutfiles[key],
            dirlist = 'lists_2017UL/'+samples[run],
            tag = tags[JERC]
        )
        command += ' '+configs[key]

os.system(command)
