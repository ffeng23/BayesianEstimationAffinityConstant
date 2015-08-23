using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using AdaptiveRejectionSampling;
using BayesianEstimateLib;
using NelderMeadMethod;

namespace Models
{
    /// <summary>
    /// this is the template for all the models that will be optimized/fitted to estimate paramters.
    /// It contains the necessary memebers including, paramters and equations, get "Y" value based on the
    /// paramters, the functions to update the parameters, the updateFunctionDistribution function updating
    /// parameters and returning the logDistributionDelegate for the specific paramter, also need function
    /// to specify which parameter we want to hold constant. The initial values of the paramters, the prior 
    /// distribution. 
    /// We also assume the prior for all regular paramters to be flat normal/diffuse prior/noninformative, 
    ///     p(theta)~N(u,sigma)  theta is the paramter, u is the prior mean chosen based on the known informaton, most time is zero
    ///                         sigma is a big value like 1E6 in order to be flat.
    /// for parameter of variance is following a gamma distribution.
    /// </summary>
    public abstract class Model
    {
        /*public Model()
        {
            this.C_ParameterList = null;
            this.C_ListFixedParameters = null;
            this.C_X = null;
            this.C_Y = null;
            this.C_ListIndexToParamsToBeUpdate = null;
            this.C_Prior = null;
            this.C_ParameterBounds = null;
        }*/
        /// <summary>
        /// the constructor to build the model with default or initial parameter values. Most of the time, in this 
        /// class, we don't have the real values of paramters, but instead we want to fit the data and estimate the
        /// values of paramters. here, but in some case, we might need to simulate the Y values then we need to specify
        /// the parameter. but not always
        /// </summary>
        /// <param name="_param">parameter values</param>
        public Model(List<double> _param)
        {
            this.C_ParameterList = _param;
            this.C_ListFixedParameters = new List<bool>();
            for (int i = 0; i < _param.Count; i++)
            {
                this.C_ListFixedParameters.Add(false);
            }
            this.C_X = null;
            this.C_Y = null;
            this.C_ListIndexToParamsToBeUpdate = null;
            this.C_Prior = null;
            this.C_ParameterBounds = null;
        }

        /// <summary>
        /// see above, we might not need a parameter values
        /// </summary>
        /// <param name="_param">parameter values</param>
        /// <param name="_X">dependent variables, could be multi-dimension</param>
        /// <param name="_Y">independent variable, only one dimension</param>
        public Model(List<double> _param, List<List<double>> _X, List<double> _Y)
        {
            
                this.C_ParameterList = _param;
                this.C_ListFixedParameters = new List<bool>();
                for (int i = 0; i < _param.Count; i++)
                {
                    this.C_ListFixedParameters.Add(false);
                }
                this.C_X = _X;
                this.C_Y = _Y;
                this.C_ListIndexToParamsToBeUpdate = null;
                this.C_Prior = null;
                this.C_ParameterBounds = null;
        }

        /// <summary>
        /// see above, this is the constructor only sepecify the data points
        /// </summary>
        /// <param name="_X">depenedent variables, could be multi-dimension</param>
        /// <param name="_Y">independent variable, only one dimension</param>
        public Model(List<List<double>> _X, List<double> _Y)
        {
            this.C_ParameterList = null;
            this.C_ListFixedParameters = null;
            this.C_X = _X;
            this.C_Y = _Y;
            this.C_ListIndexToParamsToBeUpdate = null;
            this.C_Prior = null;
            this.C_ParameterBounds = null;
        }

        public void SetupPrior(List<List<double>> _prior)
        {
            C_Prior = _prior;
        }
        public void SetupParameterBounds(List<List<double>> _bounds)
        {
            C_ParameterBounds = _bounds;
        }
        /// <summary>
        /// this is a very important function, the trick is that sometimes we might not want to estimate/vary all the 
        /// parameters, so we might want to skip some parameters. This is a trouble, since when Gibbs sampler calling to
        /// draw populations, it only iterate through the parameter lists we provided to the caller of Gibbs sampler,
        /// so the sampler doesn't know which ones to skip, etc. Then by calling this and set which one to be updated, 
        /// the sampler doesn't have to know which to skip, the model then knows which one to provide.
        /// 
        /// For example, suppose the model has 3 paramters, we only need to update/estimate 2 and 3, but not 1. Then we can specify
        /// a list containing only {1,2} (starting at 0). then when the sampler is build we feed the sample with two paramters in the order
        /// of 2 and 3. Then when the sampler iterate through 2 and 3, it get the index of the paramter 2 (index 0) and paramter 3 
        /// (Index of 1), then it feed the index to the model 0 for paramter 2, and the model will look up the list of {1,2} mentioned
        /// above, it know for index of 0 (for parameter 2) it returns the correct conditional distribution based on this input.
        /// !!!important thing, this is not a virtual member
        /// </summary>
        /// <param name="_func">the list of funcs/parameters to be updated in the sample drawing</param>
        public void setFunctionDelegateForUpdating(List<int> _func)
        {
            this.C_ListIndexToParamsToBeUpdate  = _func;
        }

        /// <summary>
        /// for properties, we don't need a virtual member it will be virtual anyway
        /// </summary>
        public int NumberOfParameters
        {
            get { return C_X.Count; }
        }
        public int NumberOfDataPoints
        {
            get { return C_X[1].Count; }
        }

        public List<List<double>> X
        {
            get { return X; }
        }
        public List<double> Y
        {
            get { return Y; }
        }
        public List<List<double>> ParameterBounds
        {
            get { return C_ParameterBounds ; }
        }
        public List<List<double>> ParameterPrior
        {
            get { return C_Prior; }
        }

        /// <summary>
        /// this is the one used to caluculate the Y values based on the parameter input. THIS IS USED BY THE SAMPLING Algorithm to calculate the proposal distribution
        /// it is has to be called by the logConditionalDistribution!!!!
        /// </summary>
        /// <param name="_x">x</param>
        /// <param name="_param">parameters</param>
        /// <returns></returns>
        public abstract double GetYValue(List<double> _x, List<double> _param);
        /// <summary>
        /// this is the actually mathematic equation to get y basedon on x. here x could be multi-dimension
        /// </summary>
        /// <param name="_x">list of independent values</param>
        /// <returns>dependent value</returns>
        public double GetYValue(List<double> _x)
        {
            if (C_ParameterList == null)
            {
                Console.WriteLine("********ERROR********: the parameters have not been set!!!");
                return double.NaN;
            }
            return GetYValue(_x, C_ParameterList);
        }
        /// <summary>
        /// simulate the dependent variable values based on the dependent values _x, assuming a Yi=Ytrue+N(0,var)
        /// a Gaussian nose with variance of _variance
        /// </summary>
        /// <param name="_x">to simulate Ys based on the input _x</param>
        /// <param name="_variance">the gaussian noise with a variance of _variance</param>
        public abstract List<double> SimulateYValues(List<List<double>> _x, double _variance=1);
        /// <summary>
        /// simulate the dependent variable values based on observed values _x inputed with the model, assuming a Yi=Ytrue+N(0,var)
        /// a Gaussian nose with variance of _variance
        /// </summary>
        /// <param name="_variance">the gaussian noise with a variance of _variance</param>
        public abstract List<double> SimulateYValues(double _variance=1);
       
        /// <summary>
        /// this is the most important function, 1)it based on the parameterCurrent input to update the model current parameters in
        ///     the model. so this is important, because the model with current paramter will be used by the Gibbs sampler to draw 
        ///     distribution; 2) it also return the log conditional pdf to the sampler based on the index, so the Sampler can call
        ///     to draw samples. and when Sampler calls to draw it will use the current model parameters.
        /// </summary>
        /// <param name="_params">the parameter list holding the current values for all the members that is being updated</param>
        /// <param name="index">the index for the one that is being updated at this step</param>
        /// <returns>function for log conditional pdf of the current parameter (indicated by index.</returns>
        public abstract LogDistributionFuctionDelegate  updateFunctionDistribution(List<double> _paramsCurrent,  int index);

        //each individual log conditional distribution of paramter is private function specified by the specific model

        //other helper functions
        /// <summary>
        /// this is a simple helper function to set all the parameter values in the model
        /// </summary>
        /// <param name="_list"></param>
        public void SetAllParameters(List<double> _list)
        {
            this.C_ParameterList=_list;
        }

        public void SetAllXs(List<List<double>> _x)
        {
            this.C_X=_x;
        }

        public void SetAllYs(List<double> _y)
        {
            this.C_Y=_y;
        }
        /// <summary>
        /// a simple helper function to set which ones to be fixed. this is another way compared with the above
        /// discribed "setFunctionDelegateForUpdating" function with the C_ListIndexToParamsToBeUpdate list. This one
        /// simply set up a flag in the list and then telling the log conditional pdf not to update it. 
        ///
        /// </summary>
        /// <param name="_i"></param>
        public void SetFixedVariables(int _i)
        {
            this.C_ListFixedParameters[_i] = true;

        }

        //declaration of member
        protected List<double> C_ParameterList;
        
        protected List<List<double>> C_X;//a list of list in case we are having multi-dimension of x
        protected List<double> C_Y; //only one dependent variable

        protected List<int> C_ListIndexToParamsToBeUpdate; //0, sloep;1, intercept;2, var.

        protected List<bool> C_ListFixedParameters; //true, means the parameter will not be estimated/fitted, keep it unchanged
                                                //false, will be need to optimize to estimate

        protected List<List<double>> C_Prior; //holding the prior distribution parameters. we assume a Gaussial flat prior for regular parameters, but invgamma for variance parameter.
                                              //[0], mean; [1], var
                                              //[0], alpha(shape); [1], beta (scale);
        protected List<List<double>> C_ParameterBounds;
    }
}
