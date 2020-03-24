import os

command = 'bash command_HTcondor.sh'

configs={}
cutfiles={}

samples = {
    'A' : 'Run2018A-17Sep2018-v2',
    'B' : 'Run2018B-17Sep2018-v1',
    'C' : 'Run2018C-17Sep2018-v1',
    'D' : 'Run2018D-PromptReco-v2',
    'MC': 'GJet_Pt-15To6000_RunIIAutumn18MiniAOD-102X',
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
        if run == 'MC':
            cutfiles[key] = cutfiles['_'.join(['A',JERC])]
        else:
            cutfiles[key] = 'config/cutFile_Run2018{run}_{JERC}.txt'.format(
            run=run,
            JERC=JERC
        )
        configs[key] = '{cutfile} {dirlist} {tag}'.format(
            cutfile = cutfiles[key],
            dirlist = 'lists_2018/'+samples[run],
            tag = tags[JERC]
        )
        command += ' '+configs[key]

os.system(command)
