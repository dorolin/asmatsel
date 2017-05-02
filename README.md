# asmatsel
simulate the effects of divergent selection and assortative mating on patterns of genome differentiation and linkage disequilibrium


files: main.C func.C asmatsel.H; to compile, try 
`g++ -Wall -O2 -o asmatsel main.C func.C -lm -lgsl -lgslcblas` (requires GSL)

to see options, type ./asmatsel without any  arguments

examples for running a simulation:
`./asmatsel -i map2.txt -l 1600 -u 0 -a 0.9 -o run2_a09.out`
or `./asmatsel -i map2S.txt -l 1600 -u 0 -s 0.8 -o run2S_s08.out`

input files (map2.txt and map2S.txt): 
can be generated with mkinfile.R

to make plot of results:
plotAsmatselMultSets.R
(e.g., using pre-processed output data run2_run2S.RData to make Phase3.pdf)
