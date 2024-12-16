# Integration test: stericInteract_fullySeparate

* Pairs of interacting rods are generated from a template .ffea input file.

* A FFEA simulation is run for each file.

* The results of each simulation are analysed, and are marked successful 
if the nodes on opposite rods are separated by more than two radii in the
final step.

* If each simulation is successful, the test passes.
