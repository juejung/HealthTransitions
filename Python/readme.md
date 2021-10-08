Markov Switching Matrices
===============================================================================

J.Jung
jjung@towson.edu

Date: 08/07/2020


In this project we estimate ologit models:
ologit health L_health X
using MEPS-HRS-2000-2016 data. The issue is that 

 * MEPS is observed at annual
   frequency with a rotating panel design, where every 2 years the panelists are
   swapped and 

 * HRS is a 'true' panel where a cohort is followed over multiple
   years and interviewed every 2! years.

This means that predictions of probabilities after an ologit from MEPS result
in Markov switching matrices based on annual frequencies and estimates based on
HRS result in Markov switching matrices with a 2 year horizon.

We attempt to use the algorithm developed in:
J. Chhatwal, S. Jayasuriya, and E. H. Elbasha, “Changing Cycle Lengths in
State-Transition Models: Challenges and Solutions,” Medical Decision Making,
vol. 36, no. 8, pp. 952–964, 2016.

to adjust the HRS 2 year Markov switching matrices to 1 year Markov switching
matrices so we can plot them side by side by age group. We use 5 year age
groups which results in 16 groups: 20, 25, 30, ...., 85, 90, 95 year olds.



