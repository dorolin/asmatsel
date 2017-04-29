# asmatsel
simulate the effects of divergent selection and assortative mating on patterns of genome differentiation and linkage disequilibrium


program requires GSL; 
files: main.C func.C asmatsel.H; to compile, try 
`g++ -Wall -O2 -o asmatsel main.C func.C -lm -lgsl -lgslcblas`

to see options, type ./asmatsel without any  arguments

example for running a simulation:
`./asmatsel -i map2S.txt -l 1600 -u 0 -s 0.8 -o run2S_s08.out >& run2S_s08.log &`

input files: 
can be generated with mkinfile.R
(map2.txt and map2S.txt)

make pdf of results:
plotAnmatselMultSets.R
(e.g., using pre-processed output data in run2_run2S.RData)
