Here in this README file, we describe what each file is about

FitController and Model are two abstract/template classes defining things.

FitController is for running ELISA related Gibbs Sampler to work and optimize the parameter. The points is that with this FitController. it is easy to call by the outside caller.

Model is for setting up the models that can be used by Gibbs Sampling algorithms. 
FitController is another layer between model and outside caller. It is tested and used by ELISA, but need to be more tested by other such as SPR.

ELISTATFitController and ELISTATModel are for using ELISTAT method to estimate the parameters.
ELISTATQuadraticFitController and ELISTATQuadraticModel are for using quadratic mothod to estimate ELISA paramter.

SprModel are the model to run spr estimation using Gibbs sampler. This model is used in by BayesianEstimationAffinityConstant and bayestimateCon.
