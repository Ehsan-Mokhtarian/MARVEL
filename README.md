# MARVEL

"On data" folder:
MARVEL.m is the main function. 
This function takes as input a data matrix and the Markov boundaries and learns a graph that is Markov equivalent to the true causal graph.
demo.m shows how to call the MARVEL function.
CI_Tests.m is the function that checks the conditional independencies. 
In the demo example, data is generated from linear Gaussian models, and the presented CI_Tests function only works for this setting.

"Oracle CI tests" folder:
In this folder, CI tests are done by checking the d-separations in the true graph. 
