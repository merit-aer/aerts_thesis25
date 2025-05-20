from LDMX.Framework import ldmxcfg

Amass='1000'

p = ldmxcfg.Process('ana')
p.sequence = [ ldmxcfg.Analyzer.from_file('mycode/analyzers/apr28Ana_2e.cxx') ]
p.inputFiles = [ 'output/simulated_data/2e_'+Amass+'MeV_10kEvents_pflowOverlay.root' ]
p.histogramFile = 'output/cool_histograms/apr28Hist_2e'+Amass+'MeV.root'
#p.maxEvents=100