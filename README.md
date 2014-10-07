gecco_2014_changing_objectives
==============================

Example code for the GECCO 2014 paper by Jonathan E. Fieldsend and Richard M. Everson 
"Efficiently Identifying Pareto Solutions when Objective Values Change"
Genetic and Evolutionary Computation Conference,
GECCO'14, pages 605-612

Institutional repository:
https://ore.exeter.ac.uk/repository/handle/10871/15248

DOI:
http://dx.doi.org/10.1145/2576768.2598279


The function GECCO_single_link_example.m will carry out the four different simulations described in the paper, and single_link_guardian_iterative.m will maintain the sink links in a changing population iteratively using any of the protocols described in the paper (note, Matlab's optimised vector comparison is used here, however an OO appoarch using references between guarding and guarding solutions should be quicker for very large populations). 
