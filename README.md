# ensemblemodeling
Ensemble modeling is used to build accurate predictive models to answer a specific question by synthesizing the results of multiple models into a single score.
This is a tiny python environment for ensemble modeling of the metabolic systems.## ensemblemodeling.py
The python script requires modules including
*python 3
*numpy
*scipy
*parallel python (pp)
The __main__ is a example and tutrial of this environment.
## PCC6803_171110.txt
This file a example metabolic model of Synechocystis sp. PCC 6803 (PCC 6803) under the photoautotrophic condition. This file consists of Reactions, Metabolites, and Parameters parts.
###//Reactions
*ID
*Stoichiometry
*type of rate equation
*effector
*Initial metabolic flux
###//Metabolites
* ID
* Initial metabolite concentration
* type of metabolite "initial" metabolite is constant
###//Parameters
* Reaction ID
* Parameter
* Lower boundary
* Reference level (not used)
* Upper boundary
* type of distribution