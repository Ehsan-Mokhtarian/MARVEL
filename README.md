# MARVEL
# Paper: A Recursive Markov Boundary-Based Approach to Causal Structure Learning

<b>"On data" folder: </b><br>
MARVEL.m is the main function.
This function takes as input a data matrix and the Markov boundaries and learns a graph that is Markov equivalent to the true causal graph. <br>
demo.m shows how to call the MARVEL function. <br>
CI_Tests.m is the function that checks the conditional independencies. <br>
In the demo example, data is generated from linear Gaussian models, and the presented CI_Tests function only works for this setting. <br>
<br>
<b>"Oracle CI tests" folder:</b> <br>
In this folder, CI tests are done by checking the d-separations in the true graph. 
