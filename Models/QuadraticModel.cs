using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Meta.Numerics.Statistics.Distributions;
using AdaptiveRejectionSampling;
using NelderMeadMethod;

namespace Models
{
    /// <summary>
    /// this is a Quadratic model with a form of y=a*x^2+b+epsilon
    /// the parameter list contains a and b and var. 0:a and 1:b and 2:var
    /// </summary>
    public class QuadraticModel:Model
    {

        
        //see the base class for description
         public QuadraticModel(List<double> _param):base(_param)
        {
             //the following are called in the base class
            /*this.C_ParameterList = _param;
            this.C_ListFixedParameters = new List<bool>();
            for (int i = 0; i < _param.Count; i++)
            {
                this.C_ListFixedParameters.Add(false);
            }
            this.C_X = null;
            this.C_Y = null;
            this.C_ListIndexToParamsToBeUpdate = null;
            this.C_Prior = null;
             */
        }

        /// <summary>
        /// see above, we might not need a parameter values
        /// </summary>
        /// <param name="_param">parameter values</param>
        /// <param name="_X">dependent variables, could be multi-dimension</param>
        /// <param name="_Y">independent variable, only one dimension</param>
        public QuadraticModel(List<double> _param, List<List<double>> _X, List<double> _Y):base(_param, _X, _Y)
        {
            //the following are called in the base class
            /*
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
            */
        }

        /// <summary>
        /// see above, this is the constructor only sepecify the data points
        /// </summary>
        /// <param name="_X">depenedent variables, could be multi-dimension</param>
        /// <param name="_Y">independent variable, only one dimension</param>
        public QuadraticModel(List<List<double>> _X, List<double> _Y):base(_X, _Y)
        {
            //the follwoing are called in the base class
            /*
            this.C_ParameterList = null;
            this.C_ListFixedParameters = null;
            this.C_X = _X;
            this.C_Y = _Y;
            this.C_ListIndexToParamsToBeUpdate = null;
            this.C_Prior = null;
             */
        }
        /*
        /// <summary>
        /// this is the actually mathematic equation to get y basedon on x. here x could be multi-dimension
        /// </summary>
        /// <param name="_x">list of independent values</param>
        /// <returns>dependent value</returns>
        public override double GetYValue(List<double> _x)
        {
            if (C_ParameterList == null)
            {
                Console.WriteLine("********ERROR********: the parameters have not been set!!!");
                return double.NaN;
            }
            return GetYValue(_x,C_ParameterList) ;
        }*/

        /// <summary>
        /// this is the one used to caluculate the Y values based on the parameter input. THIS IS USED BY THE SAMPLING Algorithm to calculate the proposal distribution
        /// it is has to be called by the logConditionalDistribution!!!!
        /// </summary>
        /// <param name="_x">x</param>
        /// <param name="_param">parameters</param>
        /// <returns></returns>
        public override double GetYValue(List<double> _x, List<double> _param)
        {
            /*if (C_ParameterList == null)
            {
                Console.WriteLine("********ERROR********: the parameters have not been set!!!");
                return double.NaN;
            }*/
            return _param[0] * _x[0]*_x[0] + _param[1];
        
        }

        /// <summary>
        /// simulate the dependent variable values based on the dependent values _x, assuming a Yi=Ytrue+N(0,var)
        /// a Gaussian nose with variance of _variance
        /// </summary>
        /// <param name="_x">to simulate Ys based on the input _x</param>
        /// <param name="_variance">the gaussian noise with a variance of _variance</param>
        public override List<double> SimulateYValues(List<List<double>> _x, double _variance = 1)
        {
            NormalDistribution nd = new NormalDistribution(0, Math.Sqrt(_variance));
            Random rng = new Random();
            List<double> temp = new List<double>();
            for (int i = 0; i < _x.Count; i++)
            {
                temp.Add(GetYValue(_x[i])+ nd.GetRandomValue(rng));
            }
            return temp;
        }
        /// <summary>
        /// simulate the dependent variable values based on observed values _x inputed with the model, assuming a Yi=Ytrue+N(0,var)
        /// a Gaussian nose with variance of _variance
        /// </summary>
        /// <param name="_variance">the gaussian noise with a variance of _variance</param>
        public override List<double> SimulateYValues(double _variance = 1)
        {
            NormalDistribution nd = new NormalDistribution(0, Math.Sqrt(_variance));
            Random rng = new Random();
            List<double> temp = new List<double>();
            for (int i = 0; i < C_X.Count; i++)
            {
                temp.Add(GetYValue(C_X[i]) + nd.GetRandomValue(rng));
            }
            return temp;
        }

        /// <summary>
        /// this is the most important function, 1)it based on the parameterCurrent input to update the model current parameters in
        ///     the model. so this is important, because the model with current paramter will be used by the Gibbs sampler to draw 
        ///     distribution; 2) it also return the log conditional pdf to the sampler based on the index, so the Sampler can call
        ///     to draw samples. and when Sampler calls to draw it will use the current model parameters.
        /// </summary>
        /// <param name="_params">the parameter list holding the current values for all the members that is being updated</param>
        /// <param name="index">the index for the one that is being updated at this step</param>
        /// <returns>function for log conditional pdf of the current parameter (indicated by index.</returns>
        public override LogDistributionFuctionDelegate updateFunctionDistribution(List<double> _paramsCurrent, int _index)
        {
            if (C_ParameterList == null)
            {
                Console.WriteLine("****ERROR***: the parameters have been set up yet");
                return null;
            }
            if (C_ParameterList.Count < _paramsCurrent.Count)
            {
                Console.WriteLine("****ERROR***: the parameter list is shorter than the input current param list, something wrong!!");
                return null;
            }
            for (int i = 0; i < _paramsCurrent.Count; i++)
            {
                switch (this.C_ListIndexToParamsToBeUpdate[i])
                {
                    case 0://parameter a
                        this.C_ParameterList[0] = _paramsCurrent[i];
                        break;
                    case 1://parameter b
                        this.C_ParameterList[1] = _paramsCurrent[i];
                        break;
                    case 2://parameter b
                        this.C_ParameterList[2] = _paramsCurrent[i];
                        break;
                                            
                    default:
                        throw new System.Exception("unknow parameter length,s omething wrong, not pointer to function distribution is returned!");
                }
            }
            
            //now return the pointer to function delegate
            switch (this.C_ListIndexToParamsToBeUpdate[_index])
            {
                case 0:

                    return LogParamACondistionalDist;

                case 1:

                    return LogParamBConditionalDist;
                case 2:

                    return LogVarConditionalDist;
                default:
                    throw new System.Exception("unknow parameter length,s omething wrong, not pointer to function distribution is returned!");
            }

        }

        //now we need to set up the log conditional distribution for each parameter
        /// <summary>
        /// the conditional distribution for the regular parameter is 
        ///     p(a|y, b, var, x)~p(a)*p(y|a,b,x,var)
        ///         p(y|a,b,x,var) is a gaussian distr with parameters with var
        ///         p(a) is a gaussian flat distr with prior
        /// </summary>
        /// <param name="_xSlope"></param>
        /// <param name="_functionNormConstant"></param>
        /// <returns></returns>
        public double LogParamACondistionalDist(double _xA, double _functionNormConstant = 0)
        {
            double temp = 0;
            //prepare the parameter
            List<double> pm = new List<double>();
            pm.Add(_xA);
            pm.Add(this.C_ParameterList[1]);
            pm.Add(this.C_ParameterList[2]);

            for (int i = 0; i < C_X.Count; i++)
            {
                temp += (C_Y[i] - GetYValue(C_X[i], pm)) * (C_Y[i] - GetYValue(C_X[i],pm));
            }
            temp /= 2 * C_ParameterList[2]; //divided by the variance param
            temp = -1 * temp - (_xA - C_Prior[0][0]) * (_xA - C_Prior[0][0]) / (2 * C_Prior[0][1]); //the below one works equally, since _xA is the current parameter value
            //temp = -1 * temp - (C_ParameterList[0] - C_Prior[0][0]) * (C_ParameterList[0] - C_Prior[0][0]) / (2 * C_Prior[0][1]);
            return temp + _functionNormConstant;
        }

        public double LogParamBConditionalDist(double _xB, double _functionNormConstant = 0)
        {
            double temp = 0;
            //prepare the parameter
            List<double> pm = new List<double>();
            pm.Add(this.C_ParameterList[0]);pm.Add(_xB);
            
            pm.Add(this.C_ParameterList[2]);

            for (int i = 0; i < C_X.Count; i++)
            {
                temp += (C_Y[i] - GetYValue(C_X[i],pm)) * (C_Y[i] - GetYValue(C_X[i],pm));
            }
            temp /= 2 * C_ParameterList[2];//divided by variance param
            temp = -1 * temp - (_xB - C_Prior[1][0]) * (_xB - C_Prior[1][0]) / (2 * C_Prior[1][1]); ;
            return temp + _functionNormConstant;
        }

        public double LogVarConditionalDist(double _xVar, double _functionNormConstant = 0)
        {
            double temp = 0;
            //prepare the parameter
            List<double> pm = new List<double>();
            
            pm.Add(this.C_ParameterList[0]);
            pm.Add(this.C_ParameterList[1]);pm.Add(_xVar);

            for (int i = 0; i < C_X.Count; i++)
            {
                temp += (C_Y[i] - GetYValue(C_X[i],pm)) * (C_Y[i] - GetYValue(C_X[i],pm));
            }
            temp /= 2 * _xVar;

            temp = -1 * temp - C_Prior[2][1] / _xVar;
            temp = temp + (C_X.Count * (-0.5) -(C_Prior[2][0]+1)) * Math.Log(_xVar);
            return temp + _functionNormConstant;
        }

        //prior
        //for prior we always have [0],mean;[1], variance;--->regular parameter
        //           [0], alpha; [1], beta; for the 
        //
             
    }//end of class
}//end of namespace
