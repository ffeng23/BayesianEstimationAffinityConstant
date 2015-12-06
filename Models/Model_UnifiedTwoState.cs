using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using Meta.Numerics.Statistics.Distributions;
using BayesianEstimateLib;
using AccessoryLib;

using NelderMeadMethod;


namespace Models
{
    /// <summary>
    /// this model is for Unified two state model to run Gibbs sampler.
    /// the member variable C_ParameterList contains the below
    ///  
    ///     the order is 0-kon_if, 1-koff_if, 2-ka_if, 3-kd_if, 4-kon_cs, 5-koff_cs, 6-ka_cs,7-kd_cs
    ///      8-Rmax, 
    ///     in this current case, we don't use offset mode. so don't worry for now.
    ///     
    /// The C_ParameterList only contains all 10 parameters, since we need the durations for the simulation.
    /// We will figure out this by looking at the input data xs.
    /// 
    /// Also, we need to be careful here, for setParameters() function, inherited from the base class,
    /// we need to pass all 14 parameters,
    /// but for updateFunctionDistribution(), we only need to pass the necessary parameters, the first 9 WITHOUT conc
    /// this parameteers are the ones to be estimated, so this list has to be identical to C_ListIndexToParamsToBeUpdate.
    /// Since this later array contains things to be updated and the former contains the current value of paramters
    /// </summary>

    class Model_UnifiedTwoState:Model
    {

        /*
         //see the base class for description
         public Model_UnifiedTwoState(List<double> _param): base(_param)
        {
            _flag_simulation = true;//means needing to run simulation again, for whatever reason.
            //empty for now
            this._SPR_TwoState_Model = new TwoStates();
        }
        */
        
        /// <summary>
        /// see above, we might not need a parameter values
        /// the member variable C_ParameterList contains the below and Again they are the estimation parameters
        ///  
        ///     the order is 0-kon_if, 1-koff_if, 2-ka_if, 3-kd_if, 4-kon_cs, 5-koff_cs, 6-ka_cs,7-kd_cs
        ///      8-Rmax, 9-variace
        ///     in this current case, we don't use offset mode. so don't worry for now.
        ///     
        /// The C_ParameterList only contains all 10 parameters, since we need the durations for the simulation.
        /// We will figure out this by looking at the input data xs.
        /// 
        /// Also, we need to be careful here, for setParameters() function, inherited from the base class,
        /// we need to pass all 14 parameters,
        /// but for updateFunctionDistribution(), we only need to pass the necessary parameters, the first 9 WITHOUT conc
        /// this parameteers are the ones to be estimated, so this list has to be identical to C_ListIndexToParamsToBeUpdate.
        /// Since this later array contains things to be updated and the former contains the current value of paramters
        /// </summary>
        /// <param name="_param">parameter values, initial values</param>
        /// <param name="_X">dependent variables, could be multi-dimension, 0: attach; 1: detach</param>
        /// <param name="_Y">independent variable, 0: attach; 1: detach</param>
        public Model_UnifiedTwoState(List<double> _param, List<List<double>> _X, List<List<double>> _Y, double _conc):base(_param)
        {
            this.C_Conc = _conc;
            _flag_simulation = true;
            this._SPR_TwoState_Model = new TwoStates();

            //now we need to do a few things 
            //the first thing is to set up expected values, 
            this.C_X = _X;
            this.C_Y = _Y[0];
            this.C_Y_Detach = _Y[1];

            //second, we need to set parameters, so far the parameter is for estimation parameters and are good
            //and this is done in the base class.
            //this.C_ParameterList = _param;

            //now we need to extract information about the duration of two phases
            ExtractSetDurationOfPhases();


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
        protected void ExtractSetDurationOfPhases()
        {
            //assuming this inputs contains 2 phases and don't have to be the input array
            //check for whether we are enlarging the maxium on.

            double duration_attach=this.C_X[0].Max();
            double duration_detach=this.C_X[1].Max();

            //compare with the input values
            this._SPR_TwoState_Model.SetDurations(new double[]{duration_attach, duration_detach, /*SimulationSPR.DefaultDeltaT*/0.5 });
        }


        /*
        /// <summary>
        /// see above, this is the constructor only sepecify the data points
        /// </summary>
        /// <param name="_X">depenedent variables, could be multi-dimension</param>
        /// <param name="_Y">independent variable, only one dimension</param>
        public Model_UnifiedTwoState(List<List<double>> _X, List<List<double>> _Y):base()
        {
            _flag_simulation = true;
            this._SPR_TwoState_Model = new TwoStates();
            //the follwoing are called in the base class
            / *
            this.C_ParameterList = null;
            this.C_ListFixedParameters = null;
            this.C_X = _X;
            this.C_Y = _Y;
            this.C_ListIndexToParamsToBeUpdate = null;
            this.C_Prior = null;
             * /
        }*/


        /// <summary>
        /// disable these ones, since we need a different one
        /// </summary>
        /// <param name="_x"></param>
        public override void SetAllXs(List<List<double>> _x)
        {
            Console.WriteLine("this one in this class \"Model_UnifiedTwoState\" is not implemented, please don't use it");
            throw new NotImplementedException("this one in this class \"Model_UnifiedTwoState\" is not implemented, please don't use it");

            base.SetAllXs(_x);
        }
        /// <summary>
        /// disable this one.
        /// </summary>
        /// <param name="_y"></param>
        public override void SetAllYs(List<double> _y)
        {
            Console.WriteLine("this one in this class \"Model_UnifiedTwoState\" is not implemented, please don't use it");
            throw new NotImplementedException("this one in this class \"Model_UnifiedTwoState\" is not implemented, please don't use it");

            base.SetAllYs(_y);
        }

        /// <summary>
        /// this is the one used to caluculate the Y values based on the parameter input. THIS IS USED BY THE SAMPLING Algorithm to calculate the proposal distribution
        /// it is has to be called by the logConditionalDistribution!!!!
        /// </summary>
        /// <param name="_x">x is list, because generally we allow the multiple dimension of x, 
        /// but here only one (set) is asked for Y alues, meaning return only on Y values. </param>
        /// <param name="_param">parameters</param>
        /// <returns></returns>
        public override double GetYValue(List<double> _x, List<double> _param)
        {
            Console.WriteLine("this one in this class \"Model_UnifiedTwoState\" is not implemented, please don't use it");
            throw new NotImplementedException("this one in this class \"Model_UnifiedTwoState\" is not implemented, please don't use it");

        }

        /// <summary>
        /// here we need to prepare the parameter in order to do the simulation.
        /// the simulation and estimation parameters are not identifical
        /// </summary>
        /// <param name="_regression_parameters"></param>
        /// <param name="_conc"></param>
        /// <param name="series_i"></param>
        /// <param name="_sim_params"></param>
        protected void PrepareParametersForSPRSimulation(ref double[] _regression_parameters, double _conc, out double[] _sim_params)
        {
            if (_regression_parameters.Count() < 10)
                throw new InvalidDataException("the parameters don't have enough element for regression and simulation");
            _sim_params = new double[17];
            //for (int i = 0; i < this._SD.NSeries; i++)
            //{
            //double conc = this._SD.SPRDataSeries[i].Concentration;
            //*****for parameters to feed into model,
            //the order is 0-kon_if, 1-koff_if, 2-ka_if, 3-kd_if, 4-kon_cs, 5-koff_cs, 6-ka_cs,7-kd_cs
            //     8-conc, 9-Rmax, 10-R0_AB,
            //     11-R0_AB_Star, 12-R0_B, 13-R0_B_Star
            //     14-attachDuration, 15-detachDuration,16-deltaT
            //double[] modelInputParams = new double[11];
            //now read the individual parameter for this series
            //for (int j = 0; j < this._NumOfIndividualParams; j++)
            //{
            _sim_params[10] = -1;//we don't know R0_AB;
            _sim_params[11] = -1;//we don't know R0_AB_Star
            _sim_params[12] = -1;//we don't know R0_B;
            _sim_params[13] = -1;//we don't know R0_B_Star
            //}
            _sim_params[0] = _regression_parameters[0];//kon_if
            _sim_params[1] = _regression_parameters[1];//koff_if
            _sim_params[2] = _regression_parameters[2];//ka_if
            _sim_params[3] = _regression_parameters[3];//kd_if
            _sim_params[4] = _regression_parameters[4];//kon_cs
            _sim_params[5] = _regression_parameters[5];//koff_cs
            _sim_params[6] = _regression_parameters[6];//ka_cs
            _sim_params[7] = _regression_parameters[7];//kd_cs

            _sim_params[8] = _conc;//conc
            _sim_params[9] = _regression_parameters[8];//Rmax


            _sim_params[14] = -1; //duration association
            _sim_params[15] = -1; //duration dissociation
            _sim_params[16] = -1;//assume we have done this. deltaT

        }

        /// <summary>
        /// this is the method intended for numerical integration functions. We know we want to do simulation as
        /// infrequent as possible, since it takes time. so want to do things in bulk meaning we want to do simulation
        /// once if possible then get all the expected Y values all at once and get the value out and pass to 
        /// the condicitional distributions. here we want the caller to tell whether the params are changed and whether 
        /// the input arrays changed. Probably every call will have a changed params except the run of logVarianceConditional
        /// but not everytime the input array is changed. so we need to be carefully. But the first run of ALL, we don need to
        /// set up new Xs for the simulation.
        /// Most importantly, the caller need to figure out whether there are new parameters and new input Xs and
        /// then set up the correct flags.
        /// </summary>
        /// <param name="_x">x is list, because generally we allow the multiple dimension of x, 
        /// but here only one (pair, attach and detach) is asked for Y alues, meaning return only on Y values. </param>
        /// <param name="_param">parameters, so far we assuming the parameters are the one </param>
        /// <returns></returns>
        protected List<List<double>> GetYValueBulk(List<List<double >> _xs, List<double> _param, bool _newParms, bool _newInputXs)
        {
            /*if (C_ParameterList == null)
            {
                Console.WriteLine("********ERROR********: the parameters have not been set!!!");
                return double.NaN;
            }*/
            //first check for the parameter has changed or not
            /*if (_param.Count != this.C_ParameterList.Count)
            {
                Console.WriteLine("********ERROR********: the parameters have not been set!!!");
                throw new ArgumentNullException("the parameter passed in is not of same dimension of the current parameter list!");
                //return double.NaN;
            }*/

            //we don't do this, we assume the user to set up the flag when we reset/change the parameters
            //just keep the things reasonable please
            /*
            for (int i = 0; i < _param.Count; i++)
            {
                //check for the new /updated paramters
                if (_param[i] - this.C_ParameterList[i] > 1E-15 || _param[i] - this.C_ParameterList[i] < -1E-15)
                {
                    _flag_simulation = true;//need to do new simulation
                    break;

                }
            }
            */
            if (_newParms)
            {
                //need to prepare parameter before calling the fit controller for running simulation
                
                this._SPR_TwoState_Model.setParameters(_param.ToArray());
                _flag_simulation = true;
            }

            if (_newInputXs)
            {
                this._SPR_TwoState_Model.addValueSToTimeArr_attach(_xs[0]);
                this._SPR_TwoState_Model.addValueSToTimeArr_detach(_xs[1]);
                _flag_simulation = true;
            }
            List<List<double>> expVals = new List<List<double>>();
            if (_flag_simulation)
            {
                
                this._SPR_TwoState_Model.run_Attach_RK();
                this._SPR_TwoState_Model.run_Detach_RK();
                
            }
            expVals.Add(  this._SPR_TwoState_Model.GetRUValuesAttach(_xs[0]));
            expVals.Add(this._SPR_TwoState_Model.GetRUValuesDetach(_xs[1]));
            //model is OD=alpha*Ab_i/(Kd+Ab_i) +beta
            return expVals;
        
        }

        /// <summary>
        /// simulate the dependent variable values based on the dependent values _x, assuming a Yi=Ytrue+N(0,var)
        /// a Gaussian nose with variance of _variance
        /// </summary>
        /// <param name="_x">to simulate Ys based on the input _x</param>
        /// <param name="_variance">the gaussian noise with a variance of _variance</param>
        public override List<double> SimulateYValues(List<List<double>> _x, double _variance = 1)
        {
            Console.WriteLine("this one in this class \"Model_UnifiedTwoState\" is not implemented, please don't use it");
            throw new NotImplementedException("this one in this class \"Model_UnifiedTwoState\" is not implemented, please don't use it");
 
            /*
            //now first check whehter _x is in the input, otherwise 


            if (!_flag_simulation)
            {
                //nothing changed, so we are good
                return;
            }
            NormalDistribution nd = new NormalDistribution(0, Math.Sqrt(_variance));
            Random rng = new Random();
            List<double> temp = new List<double>();
            for (int i = 0; i < _x.Count; i++)
            {
                temp.Add(GetYValue(_x[i])+ nd.GetRandomValue(rng));
            }
            return temp;*/
        }
        /// <summary>
        /// simulate the dependent variable values based on observed values _x inputed with the model, assuming a Yi=Ytrue+N(0,var)
        /// a Gaussian nose with variance of _variance
        /// </summary>
        /// <param name="_variance">the gaussian noise with a variance of _variance</param>
        public override List<double> SimulateYValues(double _variance = 1)
        {
            Console.WriteLine("this one in this class \"Model_UnifiedTwoState\" is not implemented, please don't use it");
            throw new NotImplementedException("this one in this class \"Model_UnifiedTwoState\" is not implemented, please don't use it");
 
            /*
            NormalDistribution nd = new NormalDistribution(0, Math.Sqrt(_variance));
            Random rng = new Random();
            List<double> temp = new List<double>();
            for (int i = 0; i < C_X.Count; i++)
            {
                temp.Add(GetYValue(C_X[i]) + nd.GetRandomValue(rng));
            }
            return temp;*/
        }

        /// <summary>
        /// this is the most important function, 1)it based on the parameterCurrent input to 
        /// update the model current parameters in
        ///     the model. so this is important, because the model with current paramter will be used by the Gibbs sampler 
        ///     to draw 
        ///     distribution; 2) it also return the log conditional pdf to the sampler based on the index, 
        ///     so the Sampler can call
        ///     to draw samples. and when Sampler calls to draw it will use the current model parameters.
        /// </summary>
        /// <param name="_params">the parameter list holding the current VALUES for all the members 
        ///                       that is being updated. This is not the index to the parameter meters
        ///                       but current VALUES of the parameters to be updated 
        ///     contains the parameter lists. 
        ///     the order is 0-kon_if, 1-koff_if, 2-ka_if, 3-kd_if, 4-kon_cs, 5-koff_cs, 6-ka_cs,7-kd_cs
        ///     8-conc, 9-Rmax, 
        ///     ----the above are the necessary ones, we only need those for now.
        ///     10-R0_AB,
        ///     11-R0_AB_Star, 12-R0_B, 13-R0_B_Star
        ///     14-attachDuration, 15-detachDuration,16-deltaT
        /// </param>
        /// <param name="index">the index for the one that is being updated at this step
        ///     _paramsCurrent contains 9 parameters, and there are
        ///     the order is 0-kon_if, 1-koff_if, 2-ka_if, 3-kd_if, 4-kon_cs, 5-koff_cs, 6-ka_cs,7-kd_cs
        ///     / * 8-conc,* / 8-Rmax, 9-variance 
        ///     
        /// </param>
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
            if (this.C_ListIndexToParamsToBeUpdate.Count > _paramsCurrent.Count)
            {
                Console.WriteLine("****ERROR***: the index to parameter list is not longer than the input current param list, something wrong!!");
                Environment.Exit(-1);
                return null;
            }

            for (int i = 0; i < this.C_ListIndexToParamsToBeUpdate.Count; i++)//_paramsCurrent HAS to be of same length as C_ListIndexToParamsToBeUpdate
            {
                int runningIndex = C_ListIndexToParamsToBeUpdate[i];
                switch (this.C_ListIndexToParamsToBeUpdate[i])  //<=== the reason that we are using ListIndextoParamsToBeUpdate is that we want to use this index to skip some elements
                {
                    case 0://parameter kon_if
                        this.C_ParameterList[0] = _paramsCurrent[runningIndex];
                        //this.C_ListFixedParameters[0] = false;
                        break;
                    case 1://parameter koff_if
                        this.C_ParameterList[1] = _paramsCurrent[runningIndex];
                        //this.C_ListFixedParameters[1] = false;
                        break;
                    case 2://parameter ka_if
                        this.C_ParameterList[2] = _paramsCurrent[runningIndex];
                        //this.C_ListFixedParameters[2] = false;
                        break;

                    case 3://parameter kd_if
                        this.C_ParameterList[3] = _paramsCurrent[runningIndex];
                        //this.C_ListFixedParameters[3] = false;
                        break;
                    case 4://parameter kon_cs
                        this.C_ParameterList[4] = _paramsCurrent[runningIndex];
                        //this.C_ListFixedParameters[4] = false;
                        break;
                    case 5://parameter koff_cs
                        this.C_ParameterList[5] = _paramsCurrent[runningIndex];
                        //this.C_ListFixedParameters[5] = false;
                        break;
                    case 6://parameter ka_cs
                        this.C_ParameterList[6] = _paramsCurrent[runningIndex];
                        //this.C_ListFixedParameters[6] = false;
                        break;

                    case 7://parameter kd_cs
                        this.C_ParameterList[7] = _paramsCurrent[runningIndex];
                        //this.C_ListFixedParameters[7] = false;
                        break;
                    case 8://parameter Rmax
                        this.C_ParameterList[8] = _paramsCurrent[runningIndex];
                        //this.C_ListFixedParameters[8] = false;
                        break;
                    case 9://parameter variance
                        this.C_ParameterList[9] = _paramsCurrent[runningIndex];
                        //this.C_ListFixedParameters[9] = false;
                        break;
                    
                    default:
                        throw new System.Exception("unknow parameter length,s omething wrong, not pointer to function distribution is returned!");
                }
            }
            
            //now return the pointer to function delegate
            switch (/*this.C_ListIndexToParamsToBeUpdate[*/_index/*]*/)
            {
                case 0:
                    return LogKonIfCondistionalDist;

                case 1:
                    return LogKoffIfConditionalDist;
                case 2:
                    return LogKaIfConditionalDist;
                case 3:
                    return LogKdIfConditionalDist;
                case 4:
                    return LogKonCsConditionalDist;

                case 5:
                    return LogKoffCsConditionalDist;
                case 6:
                    return LogKaCsConditionalDist;
                case 7:
                    return LogKdCsConditionalDist;
                case 8:
                    return LogRmaxConditionalDist;
                case 9:
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
        public double LogKonIfCondistionalDist(double _Kon_IF, double _functionNormConstant = 0)
        {

            double temp = 0;
            //prepare the parameter
            double[] pm = this.C_ParameterList.ToArray() ;
            pm[0]=_Kon_IF;
            
            //now we need to do the simulation with the updated parameters, but assuming the input x is the same
            double[] simulateParams;
            //double [] estimationParams=pm.ToArray();
            PrepareParametersForSPRSimulation(ref pm, this.C_Conc, out simulateParams);

            //now re do the simulation
            List<List<double>> expYs=this.GetYValueBulk(this.C_X, new List<double>(simulateParams), true, false);
            //association
            for (int i = 0; i < C_X[0].Count; i++)
            {
                temp += (C_Y[i] -expYs[0][i]) * (C_Y[i] - expYs[0][i]);
            }

            //dissociation
            for (int i = 0; i < C_X[1].Count; i++)
            {
                temp += (C_Y_Detach[i] - expYs[1][i]) * (C_Y_Detach[i] - expYs[1][i]);
            }
            temp /= 2 * pm[9]; //divided by the variance param
            temp = -1 * temp - (_Kon_IF  - C_Prior[0][0]) * (_Kon_IF - C_Prior[0][0]) / (2 * C_Prior[0][1]); //the below one works equally, since _xA is the current parameter value
            //temp = -1 * temp - (C_ParameterList[0] - C_Prior[0][0]) * (C_ParameterList[0] - C_Prior[0][0]) / (2 * C_Prior[0][1]);
            return temp + _functionNormConstant;
        }

        public double LogKoffIfConditionalDist(double _Koff_IF, double _functionNormConstant = 0)
        {
            double temp = 0;
            int indexOfParam = 1;
            double[] pm = this.C_ParameterList.ToArray();
            pm[indexOfParam] = _Koff_IF;

            //now we need to do the simulation with the updated parameters, but assuming the input x is the same
            double[] simulateParams;
            //double [] estimationParams=pm.ToArray();
            PrepareParametersForSPRSimulation(ref pm, this.C_Conc, out simulateParams);

            //now re do the simulation
            List<List<double>> expYs = this.GetYValueBulk(this.C_X, new List<double>(simulateParams), true, false);
            //association
            for (int i = 0; i < C_X[0].Count; i++)
            {
                temp += (C_Y[i] - expYs[0][i]) * (C_Y[i] - expYs[0][i]);
            }

            //dissociation
            for (int i = 0; i < C_X[1].Count; i++)
            {
                temp += (C_Y_Detach[i] - expYs[1][i]) * (C_Y_Detach[i] - expYs[1][i]);
            }
            temp /= 2 * pm[9];//divided by variance param
            temp = -1 * temp - (_Koff_IF - C_Prior[indexOfParam][0]) * (_Koff_IF - C_Prior[indexOfParam][0]) / (2 * C_Prior[indexOfParam][1]); ;
            return temp + _functionNormConstant;
        }
        public double LogKaIfConditionalDist(double _Ka_IF, double _functionNormConstant = 0)
        {
            double temp = 0;
            int indexOfParam = 2;
            double[] pm = this.C_ParameterList.ToArray();
            pm[indexOfParam] = _Ka_IF;

            //now we need to do the simulation with the updated parameters, but assuming the input x is the same
            double[] simulateParams;
            //double [] estimationParams=pm.ToArray();
            PrepareParametersForSPRSimulation(ref pm, this.C_Conc, out simulateParams);

            //now re do the simulation
            List<List<double>> expYs = this.GetYValueBulk(this.C_X, new List<double>(simulateParams), true, false);
            //association
            for (int i = 0; i < C_X[0].Count; i++)
            {
                temp += (C_Y[i] - expYs[0][i]) * (C_Y[i] - expYs[0][i]);
            }

            //dissociation
            for (int i = 0; i < C_X[1].Count; i++)
            {
                temp += (C_Y_Detach[i] - expYs[1][i]) * (C_Y_Detach[i] - expYs[1][i]);
            }
            temp /= 2 * C_ParameterList[9];//divided by variance param
            temp = -1 * temp - (_Ka_IF - C_Prior[indexOfParam][0]) * (_Ka_IF - C_Prior[indexOfParam][0]) / (2 * C_Prior[indexOfParam][1]); ;
            return temp + _functionNormConstant;
        }

        public double LogKdIfConditionalDist(double _Kd_IF, double _functionNormConstant = 0)
        {
            double temp = 0;
            int indexOfParam = 3;
            double[] pm = this.C_ParameterList.ToArray();
            pm[indexOfParam] = _Kd_IF;

            //now we need to do the simulation with the updated parameters, but assuming the input x is the same
            double[] simulateParams;
            //double [] estimationParams=pm.ToArray();
            PrepareParametersForSPRSimulation(ref pm, this.C_Conc, out simulateParams);

            //now re do the simulation
            List<List<double>> expYs = this.GetYValueBulk(this.C_X, new List<double>(simulateParams), true, false);
            //association
            for (int i = 0; i < C_X[0].Count; i++)
            {
                temp += (C_Y[i] - expYs[0][i]) * (C_Y[i] - expYs[0][i]);
            }

            //dissociation
            for (int i = 0; i < C_X[1].Count; i++)
            {
                temp += (C_Y_Detach[i] - expYs[1][i]) * (C_Y_Detach[i] - expYs[1][i]);
            }
            temp /= 2 * C_ParameterList[9];//divided by variance param
            temp = -1 * temp - (_Kd_IF - C_Prior[indexOfParam][0]) * (_Kd_IF - C_Prior[indexOfParam][0]) / (2 * C_Prior[indexOfParam][1]); ;
            return temp + _functionNormConstant;
        }

        public double LogKonCsConditionalDist(double _Kon_CS, double _functionNormConstant = 0)
        {
            double temp = 0;
            int indexOfParam = 4;
            double[] pm = this.C_ParameterList.ToArray();
            pm[indexOfParam] = _Kon_CS;

            //now we need to do the simulation with the updated parameters, but assuming the input x is the same
            double[] simulateParams;
            //double [] estimationParams=pm.ToArray();
            PrepareParametersForSPRSimulation(ref pm, this.C_Conc, out simulateParams);

            //now re do the simulation
            List<List<double>> expYs = this.GetYValueBulk(this.C_X, new List<double>(simulateParams), true, false);
            //association
            for (int i = 0; i < C_X[0].Count; i++)
            {
                temp += (C_Y[i] - expYs[0][i]) * (C_Y[i] - expYs[0][i]);
            }

            //dissociation
            for (int i = 0; i < C_X[1].Count; i++)
            {
                temp += (C_Y_Detach[i] - expYs[1][i]) * (C_Y_Detach[i] - expYs[1][i]);
            }
            temp /= 2 * C_ParameterList[9];//divided by variance param
            temp = -1 * temp - (_Kon_CS - C_Prior[indexOfParam][0]) * (_Kon_CS - C_Prior[indexOfParam][0]) / (2 * C_Prior[indexOfParam][1]); ;
            return temp + _functionNormConstant;
        }

        public double LogKoffCsConditionalDist(double _Koff_CS, double _functionNormConstant = 0)
        {
            double temp = 0;
            int indexOfParam = 5;
            double[] pm = this.C_ParameterList.ToArray();
            pm[indexOfParam] = _Koff_CS;

            //now we need to do the simulation with the updated parameters, but assuming the input x is the same
            double[] simulateParams;
            //double [] estimationParams=pm.ToArray();
            PrepareParametersForSPRSimulation(ref pm, this.C_Conc, out simulateParams);

            //now re do the simulation
            List<List<double>> expYs = this.GetYValueBulk(this.C_X, new List<double>(simulateParams), true, false);
            //association
            for (int i = 0; i < C_X[0].Count; i++)
            {
                temp += (C_Y[i] - expYs[0][i]) * (C_Y[i] - expYs[0][i]);
            }

            //dissociation
            for (int i = 0; i < C_X[1].Count; i++)
            {
                temp += (C_Y_Detach[i] - expYs[1][i]) * (C_Y_Detach[i] - expYs[1][i]);
            }
            temp /= 2 * C_ParameterList[9];//divided by variance param
            temp = -1 * temp - (_Koff_CS - C_Prior[indexOfParam][0]) * (_Koff_CS - C_Prior[indexOfParam][0]) / (2 * C_Prior[indexOfParam][1]); ;
            return temp + _functionNormConstant;
        }

        public double LogKaCsConditionalDist(double _Ka_CS, double _functionNormConstant = 0)
        {
            double temp = 0;
            int indexOfParam = 6;
            double[] pm = this.C_ParameterList.ToArray();
            pm[indexOfParam] = _Ka_CS;

            //now we need to do the simulation with the updated parameters, but assuming the input x is the same
            double[] simulateParams;
            //double [] estimationParams=pm.ToArray();
            PrepareParametersForSPRSimulation(ref pm, this.C_Conc, out simulateParams);

            //now re do the simulation
            List<List<double>> expYs = this.GetYValueBulk(this.C_X, new List<double>(simulateParams), true, false);
            //association
            for (int i = 0; i < C_X[0].Count; i++)
            {
                temp += (C_Y[i] - expYs[0][i]) * (C_Y[i] - expYs[0][i]);
            }

            //dissociation
            for (int i = 0; i < C_X[1].Count; i++)
            {
                temp += (C_Y_Detach[i] - expYs[1][i]) * (C_Y_Detach[i] - expYs[1][i]);
            }
            temp /= 2 * C_ParameterList[9];//divided by variance param
            temp = -1 * temp - (_Ka_CS - C_Prior[indexOfParam][0]) * (_Ka_CS - C_Prior[indexOfParam][0]) / (2 * C_Prior[indexOfParam][1]); ;
            return temp + _functionNormConstant;
        }

        public double LogKdCsConditionalDist(double _Kd_CS, double _functionNormConstant = 0)
        {
            double temp = 0;
            int indexOfParam = 7;
            double[] pm = this.C_ParameterList.ToArray();
            pm[indexOfParam] = _Kd_CS;

            //now we need to do the simulation with the updated parameters, but assuming the input x is the same
            double[] simulateParams;
            //double [] estimationParams=pm.ToArray();
            PrepareParametersForSPRSimulation(ref pm, this.C_Conc, out simulateParams);

            //now re do the simulation
            List<List<double>> expYs = this.GetYValueBulk(this.C_X, new List<double>(simulateParams), true, false);
            //association
            for (int i = 0; i < C_X[0].Count; i++)
            {
                temp += (C_Y[i] - expYs[0][i]) * (C_Y[i] - expYs[0][i]);
            }

            //dissociation
            for (int i = 0; i < C_X[1].Count; i++)
            {
                temp += (C_Y_Detach[i] - expYs[1][i]) * (C_Y_Detach[i] - expYs[1][i]);
            }
            temp /= 2 * C_ParameterList[9];//divided by variance param
            temp = -1 * temp - (_Kd_CS - C_Prior[indexOfParam][0]) * (_Kd_CS - C_Prior[indexOfParam][0]) / (2 * C_Prior[indexOfParam][1]); ;
            return temp + _functionNormConstant;
        }
        public double LogRmaxConditionalDist(double _Rmax, double _functionNormConstant = 0)
        {
            double temp = 0;
            int indexOfParam = 8;
            double[] pm = this.C_ParameterList.ToArray();
            pm[indexOfParam] = _Rmax;

            //now we need to do the simulation with the updated parameters, but assuming the input x is the same
            double[] simulateParams;
            //double [] estimationParams=pm.ToArray();
            PrepareParametersForSPRSimulation(ref pm, this.C_Conc, out simulateParams);

            //now re do the simulation
            List<List<double>> expYs = this.GetYValueBulk(this.C_X, new List<double>(simulateParams), true, false);
            //association
            for (int i = 0; i < C_X[0].Count; i++)
            {
                temp += (C_Y[i] - expYs[0][i]) * (C_Y[i] - expYs[0][i]);
            }

            //dissociation
            for (int i = 0; i < C_X[1].Count; i++)
            {
                temp += (C_Y_Detach[i] - expYs[1][i]) * (C_Y_Detach[i] - expYs[1][i]);
            }
            temp /= 2 * C_ParameterList[9];//divided by variance param
            temp = -1 * temp - (_Rmax - C_Prior[indexOfParam][0]) * (_Rmax - C_Prior[indexOfParam][0]) / (2 * C_Prior[indexOfParam][1]); ;
            return temp + _functionNormConstant;
        }
        public double LogVarConditionalDist(double _xVar, double _functionNormConstant = 0)
        {
            double temp = 0;
            int indexOfParam = 9;
            double[] pm = this.C_ParameterList.ToArray();
            pm[indexOfParam] = _xVar;

            //now we need to do the simulation with the updated parameters, but assuming the input x is the same
            double[] simulateParams;
            //double [] estimationParams=pm.ToArray();
            PrepareParametersForSPRSimulation(ref pm, this.C_Conc, out simulateParams);

            //now re do the simulation
            List<List<double>> expYs = this.GetYValueBulk(this.C_X, new List<double>(simulateParams), true, false);
            //association
            for (int i = 0; i < C_X[0].Count; i++)
            {
                temp += (C_Y[i] - expYs[0][i]) * (C_Y[i] - expYs[0][i]);
            }

            //dissociation
            for (int i = 0; i < C_X[1].Count; i++)
            {
                temp += (C_Y_Detach[i] - expYs[1][i]) * (C_Y_Detach[i] - expYs[1][i]);
            }
            temp /= 2 * _xVar;

            temp = -1 * temp - C_Prior[indexOfParam][1] / _xVar;
            temp = temp + ((C_X[0].Count + C_X[1].Count) * (-0.5) - (C_Prior[indexOfParam][0] + 1)) * Math.Log(_xVar);
            return temp + _functionNormConstant;
        }

        //new members
        protected bool _flag_simulation;
        protected TwoStates _SPR_TwoState_Model;

        //in this case, we do attach and detach seperately,
        //the original X[0] is for the attaching phase values
        //the original X[1] is for the attaching phase values 
        //the original Y is for the attaching phase values 
        
        //detaching phase Y
        List<double> C_Y_Detach;

        double C_Conc;
    }//end of class
}
