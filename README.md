Markov Health Transitions
===============================================================================

- Author: Juergen Jung
- Towson University, Department of Economics
- Email: jjung@towson.edu
- Web: https://juejung.github.io/

This repository contains some of the replication codes for my paper:

Jung, Juergen (2021), "Estimating Transition Probabilities Between Health States Using U.S.
Longitudinal Survey Data," Empirical Economics, forthcoming.

Please cite this paper if these codes turn out to be useful for your research.

Description
===========

In this project I estimate ologit (and oprobit) models of health state
transition probabilities using a combined MEPS+HRS dataset so that the entire
lifecycle can be tracked.

The estimates are based on predictions from ologit (oprobit) models and the
Long and Freese (2014) SPost13 Stata commands are used heavily.

References:

- S. J. Long and J. Freese, Regression Models for Categorical Dependent Variables Using Stata, 3 ed. College Station, TX: Stata Press, 2014.
- Web link: https://jslsoc.sitehost.iu.edu/spost13.htm

One of the issues that arises is: 

 * MEPS is observed at annual
   frequency with a rotating panel design, where every 2 years the panelists are
   swapped and 

 * HRS is a 'true' panel where a cohort is followed over multiple
   years and interviewed every 2! years.

This means that predictions of probabilities after an ologit from MEPS result
in switching matrices based on annual frequencies and estimates based on
HRS result in switching matrices with a 2 year horizon.

I use the algorithm developed in:
J. Chhatwal, S. Jayasuriya, and E. H. Elbasha, “Changing Cycle Lengths in
State-Transition Models: Challenges and Solutions,” Medical Decision Making,
vol. 36, no. 8, pp. 952–964, 2016.

to adjust the HRS-2 year probability switching matrices to 1 year probability switching
matrices in order to make the results consistent between the two surveys.
This stochastic root method requires to impose the Markov property on the
health transition probabilities.


