from LDMX.Framework import ldmxcfg

p = ldmxcfg.Process('ana')
p.sequence = [ ldmxcfg.Analyzer.from_file('mycode/analyzers/may6Ana_1e.cxx') ]
p.inputFiles = [ 'output/may6pf_10000e8GeV.root' ]
p.histogramFile = 'output/cool_histograms/may6Hist_1e.root'
#p.maxEvents=100