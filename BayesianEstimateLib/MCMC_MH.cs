using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using Meta.Numerics.Statistics.Distributions;
using System.IO;


namespace BayesianEstimateLib
{
    public sealed  class  MCMC_MH:MCMC
    {
        /// <summary>
        /// empty contructor
        /// </summary>
        public MCMC_MH():base()
        {
            //empty for now
        }

        
        /// <summary>
        /// to do one step of MCMC MH algorithm
        /// </summary>
        protected override void MCMCStep(double prevLogLLD, int steps)
        {
            double cur_ka = _ka;
            double cur_kd = _kd;
            double cur_kM = _kM;
            double cur_conc = _conc;
            double cur_Rmax = _Rmax;
            double cur_sigma = _sigma;
            double cur_R0 = this.MCMC_r0;
            //double cur_C0 = _C0_detach;
            //double cur_D0 = _D0_detach;
            
            //don't update, _deltaT and _duration

            /*double curLamda = lamda;
            double curQ = q;
            double curLamdaEffective = lamdaEffective;
            double curQEffective = qEffective;
            double curV0 = v0;
            */
            double sd_ka = 0.05;
            double sd_kd = 0.05;//keep it constant.
            double sd_kM = 0.00;
            double sd_conc = 0.05;
            double sd_Rmax = 0.05;//keep this unchange
            double sd_sigma = 0.05;//keep this unchange
            double sd_R0 = 0.05;
            //double sd_C0 = 0.0;
            //double sd_D0 = 0.0;

            //first calculate the logLikelihood of the current parameters
            /* List<double> firstGenDivTime;
            List<double> firstGenDeathTime;
            List<double> subSequentGenDivTime;
            List<double> subSequentGenDeathTime;
            List<int> modelPredictedTotalLiveCell;
            List<int> modelPredictedTotalDeadCell;
             */

            List<double> sim_ru=new List<double> ();
            double cur_loglld;
            cur_loglld = prevLogLLD;

            if(steps==0) //first step, we need calculate this.
            {
                MC_nid.setParameters(cur_ka, cur_kd, cur_kM, cur_conc, cur_Rmax,cur_R0 );
                cur_loglld = 0;
                MC_nid.run_Attach();
                sim_ru = MC_nid.RU_Attach;
                cur_loglld = logLikelihood(MC_ru_attach, sim_ru, cur_sigma);
                MC_nid.run_Detach();
                sim_ru = MC_nid.RU_Detach;
                cur_loglld += logLikelihood(MC_ru_detach, sim_ru, cur_sigma);   
            }
            
            //update the parameters
            double next_ka = Math.Exp(Math.Log(cur_ka) + sd_ka * zRand.GetRandomValue(rng));
            double next_kd = Math.Exp(Math.Log(cur_kd) + sd_kd * zRand.GetRandomValue(rng));
            double next_kM = Math.Exp(Math.Log(cur_kM) + sd_kM * zRand.GetRandomValue(rng));
            double next_conc = Math.Exp(Math.Log(cur_conc) + sd_conc * zRand.GetRandomValue(rng));
            double next_Rmax = Math.Exp(Math.Log(cur_Rmax) + sd_Rmax * zRand.GetRandomValue(rng));
            double next_R0 = Math.Exp(Math.Log(cur_R0) + sd_R0 * zRand.GetRandomValue(rng));
            //double next_C0 = Math.Exp(Math.Log(cur_C0) + sd_C0 * zRand.GetRandomValue(rng));
            //double next_D0 = Math.Exp(Math.Log(cur_D0) + sd_D0 * zRand.GetRandomValue(rng));
            /*double nextQ = uRand2.GetRandomValue(rng3);
            double nextLamdaEffective = Math.Exp(Math.Log(curLamdaEffective) + sdLE * zRand.GetRandomValue(rng));
            double nextQEffective = uRand2.GetRandomValue(rng3);
            double nextV0 = Math.Exp(Math.Log(curV0) + sdV0 * zRand.GetRandomValue(rng));
            */
            double next_sigma = Math.Exp(Math.Log(cur_sigma) + sd_sigma * zRand.GetRandomValue(rng));
            //double nextSigmaDead = Math.Exp(Math.Log(curSigmaDead) + sdSDead * zRand.GetRandomValue(rng));
            MC_nid.setParameters(next_ka, next_kd, next_kM, next_conc, next_Rmax, next_R0);
            double next_loglld=0;
            MC_nid.run_Attach();
            sim_ru = MC_nid.RU_Attach;
            next_loglld = logLikelihood(MC_ru_attach, sim_ru, next_sigma);
            MC_nid.run_Detach( );
            sim_ru = MC_nid.RU_Detach;
            next_loglld += logLikelihood(MC_ru_detach, sim_ru, next_sigma);

            bool accept;
            if (next_loglld > cur_loglld)
            {
                accept = true;
            }
            else
            {
                double u = uRand.GetRandomValue(rng2);
                //cout<<"\tnot accepted"<<endl;
                //cout<<"\tu is "<<u<<";logu is "<<log(u)<<endl;
                if (Math.Log(u) < next_loglld - cur_loglld)
                {
                    //cout<<"\tsecond accepted"<<endl;
                    accept = true;
                }
                else
                {
                    //cout<<"\tsecond Not"<<endl;
                    accept = false;
                }
            }

            //write down posterior distribution
            //for accept we need 
            if (accept)
            {
                _kaArr.Add(next_ka);
                _kdArr.Add(next_kd);
                _kMArr.Add(next_kM);
                _concArr.Add(next_conc);
                _RmaxArr.Add(next_Rmax);
                _sigmaArr.Add(next_sigma);
                MCMC_r0Arr.Add(next_R0);
                _ka = next_ka;
                _kd = next_kd;
                _kM = next_kM;
                _conc = next_conc;
                _Rmax = next_Rmax;
                _sigma = next_sigma;
                MCMC_r0 = next_R0;
                _cur_LogLld = next_loglld ;
                
            }
            else
            {
                _kaArr.Add(cur_ka);
                _kdArr.Add(cur_kd);
                _kMArr.Add(cur_kM);
                _concArr.Add(cur_conc);
                _RmaxArr.Add(cur_Rmax);
                _sigmaArr.Add(cur_sigma);
                MCMC_r0Arr.Add(cur_R0);
                _cur_LogLld = cur_loglld;
            }
            _curLldArr.Add(cur_loglld);
            _nextLldArr.Add(next_loglld);

            //save next paramters too
            writerNextParameter.WriteLine(steps +"\t"+next_ka  + "\t" + next_kd + "\t" +next_kM +"\t"
                    + next_conc +"\t"+next_Rmax  +"\t"+ next_sigma +"\t"+next_R0+"\t"+cur_loglld +"\t"+ next_loglld );

            if(steps==MC_total_cycles-1)
            {
                nextParameters.Add("conc", next_conc);
                nextParameters.Add("ka", next_ka);
                nextParameters.Add("kd", next_kd);
                nextParameters.Add("kM", next_kM);
                nextParameters.Add("Rmax", next_Rmax);
                nextParameters.Add("sigma", next_sigma);
                nextParameters.Add("R0", next_R0);

             }
        }//end of MCMCStep()

    }//end of class
}
