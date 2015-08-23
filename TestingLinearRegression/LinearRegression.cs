using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Meta.Numerics.Statistics.Distributions;
using AdaptiveRejectionSampling;
using NelderMeadMethod;

namespace TestingLinearRegression
{
    /// <summary>
    /// this is the class specify the mathmatical model for the equation we want to estimate the parameters based
    /// on the data observed. Here, we need to specify paramters, the function to generate the "Y" value basedon the 
    /// paramters and dependents. We also need to specify the updating function that could update the paramters basedon the
    /// input, also return the log of Pdf for the specific paramters.
    /// </summary>
    public class LinearRegression
    {
        public LinearRegression(double _slope, double _intercept, double _var)
        {
            C_Slope = _slope;
            C_Intercept = _intercept;
            C_Var = _var;
            C_Y = new List<double>();
           
        }
        public void setFunctionDelegateForUpdating(List<int> _func)
        {
            this.C_ListIndexToParamsToBeUpdate  = _func;
        }

        public void GetYValues(List<double> _x)
        {
            this.C_X = _x;
            NormalDistribution nd=new NormalDistribution(0, Math.Sqrt(C_Var));
            Random rng = new Random(AccessoryLib.AceessoryLib.SEED);
            for(int i=0;i<_x.Count;i++)
            {
                C_Y.Add(C_Slope * C_X[i] + C_Intercept + nd.GetRandomValue(rng));
            }
        }
        
        
        /// <summary>
        /// the order is [slope:0, intercept:1, var:2]
        /// </summary>
        /// <param name="_params">the parameter list holding the current values for all the members that is being updated</param>
        /// <param name="index">the index for the one that is being updated at this step</param>
        /// <returns></returns>
        public LogDistributionFuctionDelegate  updateFunctionDistribution(List<double> _paramsCurrent,  int index)
        {
            for (int i = 0; i < _paramsCurrent.Count; i++)
            {
                switch (this.C_ListIndexToParamsToBeUpdate[i])
                {
                    case 0:
                        C_Slope = _paramsCurrent[i];
                        break;
                    case 1:
                        C_Intercept = _paramsCurrent[i];
                        break;
                    case 2:
                        C_Var = _paramsCurrent[i];
                        break;
                    default:
                        throw new System.Exception("unknow parameter length,s omething wrong, not pointer to function distribution is returned!");
                }
            }
            /*
            C_Slope = _params[0];
            C_Intercept = _params[1];
            C_Var = _params[2];*/
            //now return the pointer to function delegate
            switch (this.C_ListIndexToParamsToBeUpdate [index])
            {
                case 0:
                   
                    return LogSlopeCondistionalDist;
                    
                case 1:
                    
                    return LogInterceptConditionalDist;
                case 2:
                    
                    return LogVarConditionalDist;
                default:
                    throw new System.Exception("unknow parameter length,s omething wrong, not pointer to function distribution is returned!");
            }
        }
        
        public double LogSlopeCondistionalDist(double _xSlope, double _functionNormConstant=0)
        {
            double temp=0;
            for(int i=0;i<C_X.Count;i++)
            {
                temp+=(C_Y[i]-(_xSlope*C_X[i]+this.C_Intercept))*(C_Y[i]-(_xSlope*C_X[i]+this.C_Intercept));
            }
            temp /=2*C_Var;
            temp = -1 * temp - (_xSlope) * (_xSlope) /( 2 * 10000);
            return temp+_functionNormConstant;
        }
        public double LogInterceptConditionalDist(double _xIntercept,double _functionNormConstant=0)
        {
            double temp = 0;
            for (int i = 0; i < C_X.Count; i++)
            {
                temp += (C_Y[i] - (C_Slope  * C_X[i] + _xIntercept)) * (C_Y[i] - (C_Slope * C_X[i] + _xIntercept));
            }
            temp /= 2 * C_Var;
            temp = -1 * temp - (_xIntercept) * (_xIntercept) / (2 * 1E10);
            return temp+_functionNormConstant;
        }

        public double LogVarConditionalDist(double _xVar, double _functionNormConstant=0)
        {
            double temp = 0;
            for (int i = 0; i < C_X.Count; i++)
            {
                temp += (C_Y[i] - (C_Slope * C_X[i] + C_Intercept)) * (C_Y[i] - (C_Slope * C_X[i] + C_Intercept));
            }
            temp /= 2 * _xVar;
            
            temp = -1 * temp - 0.001/ _xVar;
            temp = temp + (C_X.Count * (-0.5) - 1.001)*Math.Log(_xVar) ;
            return temp+ _functionNormConstant;
        }

        //declaration of member
        double C_Slope;
        double C_Intercept;
        double C_Var;

        public List<double> C_X;
        public List<double> C_Y;

        public List<int> C_ListIndexToParamsToBeUpdate; //0, sloep;1, intercept;2, var.
    }
}
