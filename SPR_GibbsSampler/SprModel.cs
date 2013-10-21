using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using AdaptiveRejectionSampling;
using BayesianEstimateLib;
using NelderMeadMethod;

namespace SPR_GibbsSampler
{
    
    public class SprModel
    {   const double deltaT=0.1;
        public SprModel(/*double _ka, double _kd, double _kM, double _conc, double _Rmax, double _R0, double _var,*/ 
            List<double> _attachT, List<double> _attachRU, List<double> _detachT, List<double> _detachRU)
        {
            /*this.C_ka = _ka;
            this.C_kd = _kd;
            this.C_kM = _kM;
            this.C_conc = _conc;
            this.C_Rmax = _Rmax;
            this.C_R0 = _R0;
            this.C_var = _var;
            */

            this.C_attachT = _attachT;
            this.C_attachRU = _attachRU;
            this.C_detachT = _detachT;
            this.C_detachRU = _detachRU;

            //now setting up the marginal distribution parameters for all the estimators

        }
        public void setFunctionDelegateForUpdating(List<int> _lst)
        {
            this.C_ListIndexToParamsToBeUpdate = _lst;
        }

        /// <summary>
        /// the order is ka:0, kd:1, kM:2, conc:3, Rmax:4, R0:5, var:6
        /// </summary>
        /// <param name="_params">this holds the current updated parameter values</param>
        /// <param name="index"></param>
        /// <returns></returns>
        public LogDistributionFuctionDelegate updateFunctionDistribution(List<double> _paramsCurrent, int index)
        {
            for (int i = 0; i < _paramsCurrent.Count; i++)
            {
                switch (this.C_ListIndexToParamsToBeUpdate[i])
                {
                    case 0:
                        C_ka = _paramsCurrent[i];
                        break;
                    case 1:
                        C_kd = _paramsCurrent[i];
                        break;
                    case 2:
                        C_kM = _paramsCurrent[i];
                        break;
                    case 3:
                        C_conc = _paramsCurrent[i];
                        break;
                    case 4:
                        C_Rmax = _paramsCurrent[i];
                        break;
                    case 5:
                        C_R0 = _paramsCurrent[i];
                        break;
                    case 6:
                        C_var = _paramsCurrent[i];
                        break;
                    default:
                        throw new System.Exception("unknow parameter length,s omething wrong, not pointer to function distribution is returned!");
                }
            }/*
            this.C_ka = _params[0];
            this.C_kd = _params[1];
            this.C_kM = _params[2];
            this.C_conc = _params[3];
            this.C_Rmax = _params[4];
            this.C_R0 = _params[5];
            this.C_var = _params[6];*/
            if (this.C_nid == null)
            {
                this.C_nid = new NumericalIntegrationOfDynamics(this.C_ka, this.C_kd, this.C_conc, this.C_Rmax,
                    this.C_R0, this.C_attachT[this.C_attachT.Count - 1] + 10, this.C_detachT[this.C_detachT.Count-1] + 10, this.C_kM,deltaT );

                //also need to update the attach/detch time array
                this.C_nid.addValueSToTimeArr_attach(this.C_attachT);
                this.C_nid.addValueSToTimeArr_detach(this.C_detachT);
            }
            
            switch (this.C_ListIndexToParamsToBeUpdate[index])
            {
                case 0:
                    this.C_updateAttach_var = true;
                    this.C_updateDetach_var = true;
                    return LogkaConditionalDist ;

                case 1:
                    this.C_updateAttach_var = true;
                    this.C_updateDetach_var = true;
                    return LogkdConditionalDist;
                case 2:
                    if (!this.C_kMFixed)
                    {
                        this.C_updateAttach_var = true;
                        this.C_updateDetach_var = true;
                    }
                    return LogkMConditionalDist;
                case 3:
                    this.C_updateAttach_var = true;
                    this.C_updateDetach_var = true;
                    return LogconcConditionalDist;
                case 4:
                    this.C_updateAttach_var = true;
                    this.C_updateDetach_var = true;
                    return LogRmaxConditionalDist;
                case 5:
                    this.C_updateAttach_var = true;
                    this.C_updateDetach_var = true;
                    return LogR0ConditionalDist;
                case 6:

                    return LogvarConditionalDist;
                default:
                    throw new System.Exception("unknow parameter length,s omething wrong, not pointer to function distribution is returned!");
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="_xka"></param>
        /// <param name="_attachAndOrdetach">_attachAndOrDetach:0, attach only; 2, detach only; 3, attach and detach</param>
        /// <param name="fixedValue">true, for fixed value that this parameter is not varying;false, varying.</param>
        /// <returns></returns>
        public double LogkaConditionalDist(double _xka, double _functionNormConstant)
        {
            if (this.C_kaFixed)
            {
                return Math.Log(this.C_ka);
            }
            
            double temp = 0;
            //now need to run numerical solution to get the y values
            C_nid.setParameters(_xka, this.C_kd, this.C_kM, this.C_conc, this.C_Rmax, this.C_R0);
            if (this.C_phaseOfEstimation == 0 || this.C_phaseOfEstimation == 3)
            {
                //this.C_nid.addValueSToTimeArr_attach(this.C_attachT);
                this.C_nid.run_Attach();
                List<double> RUValues=this.C_nid.GetRUValuesAttach(this.C_attachT);
                for (int i = 0; i < C_attachT.Count; i++)
                {
                    temp += (C_attachRU[i] - RUValues[i]) * (C_attachRU[i] - RUValues[i]);
                }
            }
            if (this.C_phaseOfEstimation == 1 || this.C_phaseOfEstimation == 3)
            {
                //this.C_nid.addValueSToTimeArr_detach(this.C_detachT);
                this.C_nid.run_Detach();
                List<double> RUValues = this.C_nid.GetRUValuesDetach(this.C_detachT);
                for (int i = 0; i < C_detachT.Count; i++)
                {
                    temp += (C_detachRU[i] - RUValues[i]) * (C_detachRU[i] - RUValues[i]);
                }
            }
            /*
            for (int i = 0; i < C_X.Count; i++)
            {
                temp += (C_Y[i] - (_xSlope * C_X[i] + this.C_Intercept)) * (C_Y[i] - (_xSlope * C_X[i] + this.C_Intercept));
            }*/
            temp /= 2 * C_var;
            temp = -1 * temp - (_xka-this.C_kaMu) * (_xka-this.C_kaMu) / (2 * this.C_kaVar );
            return temp + _functionNormConstant;
        }
        public double LogkdConditionalDist(double _xkd, double _functionNormConstant)
        {
            if (this.C_kdFixed )
            {
                return Math.Log(this.C_kd);
            }

            double temp = 0;
            //now need to run numerical solution to get the y values
            C_nid.setParameters(this.C_ka, _xkd, this.C_kM, this.C_conc, this.C_Rmax, this.C_R0);

            if (this.C_phaseOfEstimation == 0 || this.C_phaseOfEstimation == 3)
            {
                //this.C_nid.addValueSToTimeArr_attach(this.C_attachT);
                this.C_nid.run_Attach();
                List<double> RUValues = this.C_nid.GetRUValuesAttach(this.C_attachT);
                for (int i = 0; i < C_attachT.Count; i++)
                {
                    temp += (C_attachRU[i] - RUValues[i]) * (C_attachRU[i] - RUValues[i]);
                }
            }
            if (this.C_phaseOfEstimation == 1 || this.C_phaseOfEstimation == 3)
            {
                //this.C_nid.addValueSToTimeArr_detach(this.C_detachT);
                this.C_nid.run_Detach();
                List<double> RUValues = this.C_nid.GetRUValuesDetach(this.C_detachT);
                for (int i = 0; i < C_detachT.Count; i++)
                {
                    temp += (C_detachRU[i] - RUValues[i]) * (C_detachRU[i] - RUValues[i]);
                }
            }
            /*
            for (int i = 0; i < C_X.Count; i++)
            {
                temp += (C_Y[i] - (_xSlope * C_X[i] + this.C_Intercept)) * (C_Y[i] - (_xSlope * C_X[i] + this.C_Intercept));
            }*/
            temp /= 2 * C_var;
            temp = -1 * temp - (_xkd-this.C_kdMu) * (_xkd-this.C_kdMu) / (2 * this.C_kdVar);
            return temp + _functionNormConstant;
        }

        public double LogkMConditionalDist(double _xkM, double _functionNormConstant)
        {
            if (this.C_kMFixed )
            {
                return Math.Log(this.C_kM);
            }

            double temp = 0;
            //now need to run numerical solution to get the y values
            C_nid.setParameters(this.C_ka, this.C_kd , _xkM, this.C_conc, this.C_Rmax, this.C_R0);

            if (this.C_phaseOfEstimation == 0 || this.C_phaseOfEstimation == 3)
            {
                //this.C_nid.addValueSToTimeArr_attach(this.C_attachT);
                this.C_nid.run_Attach();
                List<double> RUValues = this.C_nid.GetRUValuesAttach(this.C_attachT);
                for (int i = 0; i < C_attachT.Count; i++)
                {
                    temp += (C_attachRU[i] - RUValues[i]) * (C_attachRU[i] - RUValues[i]);
                }
            }
            if (this.C_phaseOfEstimation == 1 || this.C_phaseOfEstimation == 3)
            {
                //this.C_nid.addValueSToTimeArr_detach(this.C_detachT);
                this.C_nid.run_Detach();
                List<double> RUValues = this.C_nid.GetRUValuesDetach(this.C_detachT);
                for (int i = 0; i < C_detachT.Count; i++)
                {
                    temp += (C_detachRU[i] - RUValues[i]) * (C_detachRU[i] - RUValues[i]);
                }
            }
            
            temp /= 2 * C_var;
            temp = -1 * temp - (_xkM - this.C_kMMu) * (_xkM - this.C_kMMu) / (2 * this.C_kMVar);
            return temp + _functionNormConstant;
        }
        public double LogconcConditionalDist(double _xConc, double _functionNormConstant)
        {
            if (this.C_concFixed ||this.C_phaseOfEstimation==1)
            {
                return Math.Log(this.C_conc);
            }

            double temp = 0;
            //now need to run numerical solution to get the y values
            C_nid.setParameters(this.C_ka, this.C_kd, this.C_kM, _xConc, this.C_Rmax, this.C_R0);

            if (this.C_phaseOfEstimation == 0 || this.C_phaseOfEstimation == 3)
            {
                //this.C_nid.addValueSToTimeArr_attach(this.C_attachT);
                this.C_nid.run_Attach();
                List<double> RUValues = this.C_nid.GetRUValuesAttach(this.C_attachT);
                for (int i = 0; i < C_attachT.Count; i++)
                {
                    temp += (C_attachRU[i] - RUValues[i]) * (C_attachRU[i] - RUValues[i]);
                }
            }
            /*
            if (_attachAndOrDetach == 1 || _attachAndOrDetach == 3)
            {
                this.C_nid.addValueSToTimeArr_detach(this.C_detachT);
                this.C_nid.run_Detach();
                List<double> RUValues = this.C_nid.GetRUValuesDetach(this.C_detachT);
                for (int i = 0; i < C_detachT.Count; i++)
                {
                    temp += (C_detachRU[i] - RUValues[i]) * (C_detachRU[i] - RUValues[i]);
                }
            }*/
            
            temp /= 2 * C_var;
            temp = -1 * temp - (_xConc - this.C_concMu) * (_xConc - this.C_concMu) / (2 * this.C_concVar);
            return temp + _functionNormConstant;
        }
        public double LogRmaxConditionalDist(double _xRmax, double _functionNormConstant)
        {
            if (this.C_RmaxFixed ||this.C_phaseOfEstimation==1)
            {
                return Math.Log(this.C_Rmax);
            }

            double temp = 0;
            //now need to run numerical solution to get the y values
            C_nid.setParameters(this.C_ka, this.C_kd, this.C_kM, this.C_conc, _xRmax, this.C_R0);

            if (this.C_phaseOfEstimation == 0 || this.C_phaseOfEstimation == 3)
            {
                //this.C_nid.addValueSToTimeArr_attach(this.C_attachT);
                this.C_nid.run_Attach();
                List<double> RUValues = this.C_nid.GetRUValuesAttach(this.C_attachT);
                for (int i = 0; i < C_attachT.Count; i++)
                {
                    temp += (C_attachRU[i] - RUValues[i]) * (C_attachRU[i] - RUValues[i]);
                }
            }
            /*
            if (_attachAndOrDetach == 1 || _attachAndOrDetach == 3)
            {
                this.C_nid.addValueSToTimeArr_detach(this.C_detachT);
                this.C_nid.run_Detach();
                List<double> RUValues = this.C_nid.GetRUValuesDetach(this.C_detachT);
                for (int i = 0; i < C_detachT.Count; i++)
                {
                    temp += (C_detachRU[i] - RUValues[i]) * (C_detachRU[i] - RUValues[i]);
                }
            }*/
            temp /= 2 * C_var;
            temp = -1 * temp - (_xRmax - this.C_RmaxMu) * (_xRmax - this.C_RmaxMu) / (2 * this.C_RmaxVar);
            return temp + _functionNormConstant;
        }
        public double LogR0ConditionalDist(double _xR0, double _functionNormConstant)
        {
            if (this.C_R0Fixed ||this.C_phaseOfEstimation==0)
            {
                return Math.Log(this.C_R0);
            }

            double temp = 0;
            //now need to run numerical solution to get the y values
            C_nid.setParameters(this.C_ka, this.C_kd, this.C_kM, this.C_conc, this.C_Rmax, _xR0);
            /*
            if (_attachAndOrDetach == 0 || _attachAndOrDetach == 3)
            {
                this.C_nid.addValueSToTimeArr_attach(this.C_attachT);
                this.C_nid.run_Attach();
                List<double> RUValues = this.C_nid.GetRUValuesAttach(this.C_attachT);
                for (int i = 0; i < C_attachT.Count; i++)
                {
                    temp += (C_attachRU[i] - RUValues[i]) * (C_attachRU[i] - RUValues[i]);
                }
            }
            */
            if (this.C_phaseOfEstimation == 1 || this.C_phaseOfEstimation == 3)
            {
                //this.C_nid.addValueSToTimeArr_detach(this.C_detachT);
                this.C_nid.run_Detach();
                List<double> RUValues = this.C_nid.GetRUValuesDetach(this.C_detachT);
                for (int i = 0; i < C_detachT.Count; i++)
                {
                    temp += (C_detachRU[i] - RUValues[i]) * (C_detachRU[i] - RUValues[i]);
                }
            }
            temp /= 2 * C_var;
            temp = -1 * temp - (_xR0 - this.C_R0Mu) * (_xR0 - this.C_R0Mu) / (2 * this.C_R0Var);
            return temp + _functionNormConstant;
        }

        public double LogvarConditionalDist(double _xVar, double _functionNormConstant)
        {
            //Console.WriteLine("in var conditions..........");

            if (this.C_varFixed )
            {
                return Math.Log(this.C_var);
            }

            double temp = 0;
            //now need to run numerical solution to get the y values
            C_nid.setParameters(this.C_ka, this.C_kd, this.C_kM, this.C_conc, this.C_Rmax, this.C_R0);
            int count = 0;
            List<double> RUValues;
            if (this.C_phaseOfEstimation==0 || this.C_phaseOfEstimation==3)
            {
                //Console.WriteLine("in var phase 0...");
                //this.C_nid.addValueSToTimeArr_attach(this.C_attachT);
                if (this.C_updateAttach_var)
                {
                    //Console.WriteLine("in var updating............");
                    this.C_nid.run_Attach();
                    this.C_updateAttach_var = false;
                }
                RUValues = this.C_nid.GetRUValuesAttach(this.C_attachT);
                for (int i = 0; i < C_attachT.Count; i++)
                {
                    count++;
                    temp += (C_attachRU[i] - RUValues[i]) * (C_attachRU[i] - RUValues[i]);
                }
                
                
            }

            if (this.C_phaseOfEstimation == 1 || this.C_phaseOfEstimation == 3)
            {
                //this.C_nid.addValueSToTimeArr_detach(this.C_detachT);
                if (this.C_updateDetach_var)
                {
                    this.C_nid.run_Detach();
                    this.C_updateDetach_var = false;
                }
                RUValues = this.C_nid.GetRUValuesDetach(this.C_detachT);
                for (int i = 0; i < C_detachT.Count; i++)
                {
                    count++;
                    temp += (C_detachRU[i] - RUValues[i]) * (C_detachRU[i] - RUValues[i]);
                }
                
            }
            
            temp /= 2 * _xVar;
            temp = -1 * temp - C_varBeta  / _xVar;
            temp = temp +  (count * (-0.5) - (1.0+C_VarAlpha))*Math.Log(_xVar);

            return temp + _functionNormConstant;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="_i">0, attaching only;1, detaching only; 2, both phases.</param>
        public void SetEstimationPhaseMode(int _i)
        {
            switch(_i)
            {
                case 0:
                case 1:
                case 2:
                    this.C_phaseOfEstimation = _i;
                    break;
                default:
                    this.C_phaseOfEstimation = 0;
                    Console.WriteLine("unknow phase mode to set, using the default mode (0).");
                    break;
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="_i">0, ka;1, kd;2,kM;3,conc;4,Rmax;5, R0;6,var; please be very careful for setting "var" to be fixed, it might break the code</param>
        public void SetFixedVariables(int _i)
        {
            switch (_i)
            {
                case 0:
                    this.C_kaFixed = true;
                    break;
                case 1:
                    this.C_kdFixed = true;
                    break;
                case 2:
                    this.C_kMFixed = true;
                    break;
                case 3:
                    this.C_concFixed = true;
                    break;
                case 4:
                    this.C_RmaxFixed = true;
                    break;
                case 5:
                    this.C_R0Fixed = true;
                    break;
                case 6: //please be very careful for setting this to be fixed, it might break the code
                    this.C_varFixed = true;
                    break;
                
                default:
                    this.C_kaFixed = this.C_kdFixed = this.C_kMFixed = this.C_concFixed = this.C_RmaxFixed = this.C_R0Fixed = this.C_varFixed = false;
                    Console.WriteLine("the setting fixed variable method is not specifying the correct values for any variable to be fixed....");
                    break;
            }
        }
        public void SetAllParamters(double _ka, double _kd, double _kM, double _conc, double _Rmax, double _R0, double _var)
        {
            this.C_ka = _ka;
            this.C_kd = _kd;
            this.C_kM = _kM;
            this.C_conc = _conc;
            this.C_Rmax = _Rmax;
            this.C_R0 = _R0;
            this.C_var = _var;
        }

        //declaration of members
        double C_ka;
        double C_kd;
        double C_kM;
        double C_conc;
        double C_Rmax;
        double C_R0;
        double C_var;
        
        NumericalIntegrationOfDynamics C_nid;
        List<double> C_attachT;
        List<double> C_attachRU;
        List<double> C_detachT;
        List<double> C_detachRU;

        //marginal distribution parameters for the estimated variables 
        double C_kaMu = 0;
        double C_kaVar = 1E17;
        double C_kdMu = 0;
        double C_kdVar = 1E5;
        double C_kMMu = 0;
        double C_kMVar = 1E35;
        double C_concMu = 0;
        double C_concVar = 1E3;
        double C_RmaxMu = 0;
        double C_RmaxVar = 1E11;
        double C_R0Mu = 00;
        double C_R0Var = 1E7;
        double C_VarAlpha = 0.1;
        double C_varBeta = 0.1;

        //extra controls
        int C_phaseOfEstimation = 3;//0, for attaching only;1, detaching only; 3, both;
        bool C_kaFixed = false;
        bool C_kdFixed = false;                                                                                                                             
        bool C_kMFixed = false;
        bool C_concFixed = false;
        bool C_RmaxFixed = false;
        bool C_R0Fixed = false;
        bool C_varFixed = false;

        //update signal for the parameter affecting nid numerical solution for C_var distribution function
        bool C_updateAttach_var=true;
        bool C_updateDetach_var = true;
        List<int> C_ListIndexToParamsToBeUpdate;


    }//end of class
}//end of namespace
