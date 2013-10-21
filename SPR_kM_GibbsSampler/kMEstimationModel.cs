using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BayesianEstimateLib;
using AdaptiveRejectionSampling;

using NelderMeadMethod;

namespace SPR_kM_GibbsSampler
{
    public class kMEstimationModel
    {
        const double deltaT_Km = 0.1;
        /// <summary>
        /// the constructor 
        /// </summary>
        /// <param name="_time_arr">list of time array, could be several pieces of data to estimate together</param>
        /// <param name="_RU_arr">as above</param>
        public kMEstimationModel(List<List<double>> _attach_T, List<List<double>> _attach_RU, List<List<double>> _detach_T, List<List<double>> _detach_RU)
        {
            this.C_attachT = _attach_T;
            this.C_attachRU = _attach_RU;
            this.C_detachT = _detach_T;
            this.C_detachRU = _detach_RU;
            this.C_totalNumOfInput = _attach_RU.Count();
        }
        /// <summary>
        /// set up the kM related parameters
        /// </summary>
        /// <param name="_molecularWeight"></param>
        /// <param name="_flowRate"></param>
        public void SetkMRelatedParameters(double _molecularWeight, List<double> _flowRate, double _temp=25)
        {
            if (_flowRate.Count != C_attachRU.Count)
            {
                throw new System.Exception("the flow rates parameters are not set up correctly");
            }

            this.C_molecularWeight = _molecularWeight;
            this.C_flowRate = _flowRate;
            this.C_temp=_temp;
            this.C_D = this.CalculateDiffusionConstant();
            this.C_kM = new List<double>(this.C_flowRate.Count);
            /*for (int i = 0; i < this.C_flowRate.Count; i++)
            {
                this.C_kM.Add();
            }*/

        }
        /// <summary>
        /// set all the necessary paramters for the estimation. these parameters are used to run numerical estimation of solutions
        /// </summary>
        /// <param name="_ka"></param>
        /// <param name="_kd"></param>
        /// <param name="_kM"></param>
        /// <param name="_conc"></param>
        /// <param name="_Rmax"></param>
        /// <param name="_R0"></param>
        /// <param name="_var"></param>
        public void SetAllParamters(List<double> _ka, List<double> _kd, List<double> _conc, List<double> _Rmax, List<double> _R0, List<double> _var, double _Ckc)
        {
            this.C_ka = _ka;
            this.C_kd = _kd;
            //this.C_kM = _kM;
            this.C_conc = _conc;
            this.C_Rmax = _Rmax;
            this.C_R0 = _R0;
            this.C_var = _var;
            this.C_Ckc = _Ckc;

            //set up C_kM
            CalculatekM(_Ckc);

            //now every should be ready and we need to initialize the nid's
            
        }

        /// <summary>
        /// in this case, we only have var and Ckc to be updated!!!
        /// </summary>
        /// <param name="_lst"></param>
        public void setFunctionDelegateForUpdating(List<int> _lst)
        {
            this.C_ListIndexToParamsToBeUpdate = _lst;
        }

        /// <summary>
        /// the order is Ckc(kM):0; Rmax0:1, Rmax1:2,...,Rmax(#-1):#;
        /// R0_0:1*#+1;R0_1:1*#+2,....,R0_(#-1):1*#+#
        /// var0:2*#+1, var1:2*#+2....var(#-1):2*#+#; for example if there are two different dataset (each dataset contains both attaching and detaching data)
        /// then, we have 0 for Ckc and 1 for var0 and 2 for var1;
        /// Ckc is one constant for all the datasets, since this one is machine property related, it is constant as
        /// long as we are running exps on the same machine. so set it up as a list, but the only first element [0] is meaningful.
        /// </summary>
        /// <param name="_params">this parameter array holds the most updated parameter values updated by Gibbs sampler</param>
        /// <param name="index"></param>
        /// <returns></returns>
        public LogDistributionFuctionDelegate updateFunctionDistribution(List<double> _paramsCurrent, int _index)
        {
            for (int i = 0; i < _paramsCurrent.Count; i++)
            {//we updated all of them 
                if(this.C_ListIndexToParamsToBeUpdate[i]==0)
                {
                        C_Ckc= _paramsCurrent[i];
                        CalculatekM(C_Ckc);
                }
                else
                {
                    if (this.C_ListIndexToParamsToBeUpdate[i] <= C_totalNumOfInput)
                    {
                        //Rmax
                        C_Rmax[this.C_ListIndexToParamsToBeUpdate[i] - 1] = _paramsCurrent[i]; 
                    }
                    else
                    {
                        if (this.C_ListIndexToParamsToBeUpdate[i] <= C_totalNumOfInput*2)
                        {
                            //R0
                            C_R0[this.C_ListIndexToParamsToBeUpdate[i] - this.C_totalNumOfInput - 1] = _paramsCurrent[i];
                        }
                        else
                        {
                            //here we assume the var are different for different dataset and the index for different var is index-1
                            C_var[this.C_ListIndexToParamsToBeUpdate[i]-this.C_totalNumOfInput*2 - 1] = _paramsCurrent[i];
                        }
                    }
                        
                }
            }
            /*
            this.C_ka = _params[0];
            this.C_kd = _params[1];
            this.C_kM = _params[2];
            this.C_conc = _params[3];
            this.C_Rmax = _params[4];
            this.C_R0 = _params[5];
            this.C_var = _params[6];*/
            //this.SetAllParamters(this.C_ka, this.C_kd, this.C_conc, this.C_Rmax, this.C_R0, this.C_var, this.C_Ckc);
            if (this.C_nid == null||this.C_nid.Count==0)
            {
                this.C_nid = new List<NumericalIntegrationOfDynamics>(this.C_flowRate.Count);
                for (int i = 0; i < this.C_flowRate.Count; i++)
                {
                    C_nid.Add(new NumericalIntegrationOfDynamics(this.C_ka[i], this.C_kd[i], this.C_conc[i], this.C_Rmax[i],
                            this.C_R0[i], this.C_attachT[i][this.C_attachT[i].Count - 1] + 1, this.C_detachT[i][this.C_detachT[i].Count - 1] + 1,
                            this.C_kM[i], deltaT_Km)
                            );
                
                //also need to update the attach/detch time array
                this.C_nid[i].addValueSToTimeArr_attach(this.C_attachT[i]);
                this.C_nid[i].addValueSToTimeArr_detach(this.C_detachT[i]);
                }
            }
            this.C_CurrentIndexOfParameterBeingUpdated = _index;

           if(this.C_ListIndexToParamsToBeUpdate[_index]==0)
            {
                
                    this.C_updateAttach_var = true;
                    this.C_updateDetach_var = true;
                    return LogCkcConditionalDist;
           }
           else
           {
               if (this.C_ListIndexToParamsToBeUpdate[_index] <= this.C_totalNumOfInput)
               {
                   this.C_updateAttach_var = true;
                   this.C_updateDetach_var = true;
                   return LogRmaxConditionalDist;
               }
               else
               {
                   if (this.C_ListIndexToParamsToBeUpdate[_index] <= this.C_totalNumOfInput * 2)
                   {
                       this.C_updateAttach_var = true;
                       this.C_updateDetach_var = true;
                       return LogR0ConditionalDist;
                   }
                   else
                   {
                       this.C_updateAttach_var = true;
                       this.C_updateDetach_var = true;
                       return LogvarConditionalDist;
                   }
               }
           }
        }
        public double LogCkcConditionalDist(double _xCkc, double _functionNormConstant)
        {
            /*if (this.C_kMFixed)
            {
                return Math.Log(this.C_kM);
            }*/

            double temp = 0;
            double total_temp = 0;
            //now need to run numerical solution to get the y values
            for (int j = 0; j < this.C_nid.Count; j++)
            {
                temp = 0;
                CalculatekM(_xCkc);
                C_nid[j].setParameters(this.C_ka[j], this.C_kd[j], this.C_kM[j], this.C_conc[j], this.C_Rmax[j], this.C_R0[j]);

                if (this.C_phaseOfEstimation == 0 || this.C_phaseOfEstimation == 3)
                {
                    //this.C_nid.addValueSToTimeArr_attach(this.C_attachT);
                    this.C_nid[j].run_Attach();
                    List<double> RUValues = this.C_nid[j].GetRUValuesAttach(this.C_attachT[j]);
                    for (int i = 0; i < C_attachT[j].Count; i++)
                    {
                        temp += (C_attachRU[j][i] - RUValues[i]) * (C_attachRU[j][i] - RUValues[i]);
                    }
                }
                if (this.C_phaseOfEstimation == 1 || this.C_phaseOfEstimation == 3)
                {
                    //this.C_nid.addValueSToTimeArr_detach(this.C_detachT);
                    this.C_nid[j].run_Detach();
                    List<double> RUValues = this.C_nid[j].GetRUValuesDetach(this.C_detachT[j]);
                    for (int i = 0; i < C_detachT[j].Count; i++)
                    {
                        temp += (C_detachRU[j][i] - RUValues[i]) * (C_detachRU[j][i] - RUValues[i]);
                    }
                }
                 
                temp /= 2 * C_var[j];
                temp = -1 * temp - (_xCkc - this.C_CkcMu) * (_xCkc- this.C_CkcMu) / (2 * this.C_CkcVar);
                total_temp += temp;
            }
            return total_temp + _functionNormConstant;
        }
        public double LogRmaxConditionalDist(double _xRmax, double _functionNormConstant)
        {
            /*if (this.C_kMFixed)
            {
                return Math.Log(this.C_kM);
            }*/

            double temp = 0;
            //double total_temp = 0;
            //now need to run numerical solution to get the y values
                            temp = 0;
                //CalculatekM(this.C_Ckc);
                int currentIndex = this.C_ListIndexToParamsToBeUpdate[this.C_CurrentIndexOfParameterBeingUpdated] - 1;
                C_nid[currentIndex].setParameters(this.C_ka[currentIndex], this.C_kd[currentIndex], this.C_kM[currentIndex], this.C_conc[currentIndex], _xRmax, this.C_R0[currentIndex]);

                if (this.C_phaseOfEstimation == 0 || this.C_phaseOfEstimation == 3)
                {
                    //this.C_nid.addValueSToTimeArr_attach(this.C_attachT);
                    this.C_nid[currentIndex].run_Attach();
                    List<double> RUValues = this.C_nid[currentIndex].GetRUValuesAttach(this.C_attachT[currentIndex]);
                    for (int i = 0; i < C_attachT[currentIndex].Count; i++)
                    {
                        temp += (C_attachRU[currentIndex][i] - RUValues[i]) * (C_attachRU[currentIndex][i] - RUValues[i]);
                    }
                }
            /*
                if (this.C_phaseOfEstimation == 1 || this.C_phaseOfEstimation == 3)
                {
                    //this.C_nid.addValueSToTimeArr_detach(this.C_detachT);
                    this.C_nid[currentIndex].run_Detach();
                    List<double> RUValues = this.C_nid[currentIndex].GetRUValuesDetach(this.C_detachT[currentIndex]);
                    for (int i = 0; i < C_detachT[currentIndex].Count; i++)
                    {
                        temp += (C_detachRU[currentIndex][i] - RUValues[i]) * (C_detachRU[currentIndex][i] - RUValues[i]);
                    }
                }
            */
                temp /= 2 * C_var[currentIndex];
                temp = -1 * temp - (_xRmax - this.C_RmaxMu) * (_xRmax - this.C_RmaxMu) / (2 * this.C_RmaxVar);
                //total_temp += temp;
            //}
            return temp + _functionNormConstant;
        }

        public double LogR0ConditionalDist(double _xR0, double _functionNormConstant)
        {
            /*if (this.C_kMFixed)
            {
                return Math.Log(this.C_kM);
            }*/

            double temp = 0;
            //double total_temp = 0;
            //now need to run numerical solution to get the y values
            temp = 0;
            //CalculatekM(this.C_Ckc);
            int currentIndex = this.C_ListIndexToParamsToBeUpdate[this.C_CurrentIndexOfParameterBeingUpdated] - C_totalNumOfInput - 1;
            C_nid[currentIndex].setParameters(this.C_ka[currentIndex], this.C_kd[currentIndex], this.C_kM[currentIndex], this.C_conc[currentIndex], this.C_Rmax[currentIndex], _xR0);
            /*
            if (this.C_phaseOfEstimation == 0 || this.C_phaseOfEstimation == 3)
            {
                //this.C_nid.addValueSToTimeArr_attach(this.C_attachT);
                this.C_nid[currentIndex].run_Attach();
                List<double> RUValues = this.C_nid[currentIndex].GetRUValuesAttach(this.C_attachT[currentIndex]);
                for (int i = 0; i < C_attachT[currentIndex].Count; i++)
                {
                    temp += (C_attachRU[currentIndex][i] - RUValues[i]) * (C_attachRU[currentIndex][i] - RUValues[i]);
                }
            }*/
            if (this.C_phaseOfEstimation == 1 || this.C_phaseOfEstimation == 3)
            {
                //this.C_nid.addValueSToTimeArr_detach(this.C_detachT);
                this.C_nid[currentIndex].run_Detach();
                List<double> RUValues = this.C_nid[currentIndex].GetRUValuesDetach(this.C_detachT[currentIndex]);
                for (int i = 0; i < C_detachT[currentIndex].Count; i++)
                {
                    temp += (C_detachRU[currentIndex][i] - RUValues[i]) * (C_detachRU[currentIndex][i] - RUValues[i]);
                }
            }

            temp /= 2 * C_var[currentIndex];
            temp = -1 * temp - (_xR0 - this.C_R0Mu) * (_xR0 - this.C_R0Mu) / (2 * this.C_R0Var);
            //total_temp += temp;
            //}
            return temp + _functionNormConstant;
        }

        /// <summary>
        /// here, this implementation is different than the SPRModel one, it assume the var is different from data set to dataset, so
        /// here we need to sample the different var sepearately. this means we don't need to draw sample from other dataset if we only 
        /// do it for one sepecific dataset. check to see which one we are sampling from (C_ListIndexToParamsToBeUpdate[this.C_CurrentIndexOfParameterBeingUpdated])
        /// </summary>
        /// <param name="_xVar"></param>
        /// <param name="_functionNormConstant"></param>
        /// <returns></returns>
        public double LogvarConditionalDist(double _xVar, double _functionNormConstant)
        {
            //Console.WriteLine("in var conditions..........");

            double temp = 0;
            //double total_temp=0;
            //for (int j = 0; j < this.C_nid.Count; j++)
            //{
                temp = 0;
                //CalculatekM(this.C_Ckc);
                int currentIndex = this.C_ListIndexToParamsToBeUpdate[this.C_CurrentIndexOfParameterBeingUpdated]- C_totalNumOfInput*2 - 1;
                C_nid[currentIndex].setParameters(this.C_ka[currentIndex], this.C_kd[currentIndex], this.C_kM[currentIndex], this.C_conc[currentIndex], this.C_Rmax[currentIndex], this.C_R0[currentIndex]);

                int count = 0;
                List<double> RUValues;
                if (this.C_phaseOfEstimation == 0 || this.C_phaseOfEstimation == 3)
                {
                    //Console.WriteLine("in var phase 0...");
                    //this.C_nid.addValueSToTimeArr_attach(this.C_attachT);
                    if (this.C_updateAttach_var)
                    {
                        //Console.WriteLine("in var updating............");
                        this.C_nid[currentIndex].run_Attach();
                        this.C_updateAttach_var = true;  //here we try update it every round to make the code easier!!!!! be careful!!!
                    }
                    RUValues = this.C_nid[currentIndex].GetRUValuesAttach(this.C_attachT[currentIndex]);
                    for (int i = 0; i < C_attachT[currentIndex].Count; i++)
                    {
                        count++;
                        temp += (C_attachRU[currentIndex][i] - RUValues[i]) * (C_attachRU[currentIndex][i] - RUValues[i]);
                    }


                }

                if (this.C_phaseOfEstimation == 1 || this.C_phaseOfEstimation == 3)
                {
                    //this.C_nid.addValueSToTimeArr_detach(this.C_detachT);
                    if (this.C_updateDetach_var)
                    {
                        this.C_nid[currentIndex].run_Detach();
                        this.C_updateDetach_var = true; //here we try update it every round to make the code easier!!!!! be careful!!!
                    }
                    RUValues = this.C_nid[currentIndex].GetRUValuesDetach(this.C_detachT[currentIndex]);
                    for (int i = 0; i < C_detachT[currentIndex].Count; i++)
                    {
                        count++;
                        temp += (C_detachRU[currentIndex][i] - RUValues[i]) * (C_detachRU[currentIndex][i] - RUValues[i]);
                    }

                }

                temp /= 2 * _xVar;
                temp = -1 * temp - C_varBeta / _xVar;
                temp = temp + (count * (-0.5) - (1.0 + C_VarAlpha)) * Math.Log(_xVar);
                //total_temp += temp;
            //}
            return temp + _functionNormConstant;
        }
        private void CalculatekM(double _Ckc)
        {
            this.C_kM = new List<double>(this.C_flowRate.Count);
            
            this.C_D = CalculateDiffusionConstant();
            for (int i = 0; i < this.C_flowRate.Count; i++)
            {
                this.C_kM.Add(_Ckc*Math.Pow(this.C_D*this.C_D*this.C_flowRate[i], 1.0/3.0));
            }
        }
        private double CalculateDiffusionConstant()
        {
            double k=1.381E-23; //#boltzmann's constant
            double T = this.C_temp + 273.15; //#absolute Temp at 25oC

            double v=7.3E-4;//  #specific density normall 7.0~7.55E-4 for protein
            double Na = 6.022E23; //  # Avogardro's number
            double eta =0.001;  //#viscosity
            double flf0=1.2; //# friction factor 1.2 a reasonable guess for most globular protein
            
            double D=k*T/( 6*Math.PI*Math.Pow((3*this.C_molecularWeight*v/(4*Math.PI*Na)),(1.0/3.0))*eta*flf0);
            return D;

            /*###reference: Christensen, L.L., Theoretical analysis of protein concentration determination using biosensor technology under conditions of partial mass transport limitation. Analytical biochemistry, 1997. 249(2): p. 153-64.
            ####25oC
            #MW: 150KD D is  5.179461e-12
            #MW: 300KD D is  4.110941e-12
            #MW: 75KD D is  6.525711e-12

            ####20oC
            #MW: 150KD D is  5.092601e-12

            ####################**************/
        }

        //********************************
        //member declaration
        //declaration of members
        List<double> C_ka;
        List<double> C_kd;
        List<double> C_kM;
        List<double> C_conc;
        List<double> C_Rmax;
        List<double> C_R0;
        List<double> C_var;
        double C_Ckc;

        List<NumericalIntegrationOfDynamics> C_nid;


        //kM related parameters
        double C_molecularWeight;
        List<double> C_flowRate;
        double C_temp;
        double C_D;


        List<List<double>> C_attachT;
        List<List<double>> C_attachRU;
        List<List<double>> C_detachT;
        List<List<double>> C_detachRU;

        //update signal for the parameter affecting nid numerical solution for C_var distribution function
        bool C_updateAttach_var = true;
        bool C_updateDetach_var = true;
        List<int> C_ListIndexToParamsToBeUpdate;
        int C_phaseOfEstimation = 3;

        //marginal distribution parameters for the estimated variables 
        /*double C_kaMu = 0;
        double C_kaVar = 1E15;
        double C_kdMu = 0;
        double C_kdVar = 1E2;
        double C_kMMu = 0;
        double C_kMVar = 1E15;
        double C_concMu = 0;
        double C_concVar = 1E2;*/
        double C_RmaxMu = 0;
        double C_RmaxVar = 1E4;
        double C_R0Mu = 00;
        double C_R0Var = 1E4;
        double C_VarAlpha = 0.1;
        double C_varBeta = 0.1;
        double C_CkcMu = 0;
        double C_CkcVar = 1E30;

        int C_CurrentIndexOfParameterBeingUpdated;
        int C_totalNumOfInput;
    }
}
