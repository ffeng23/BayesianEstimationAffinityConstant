using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BayesianEstimateLib
{
    public class MCMC_Gibbs:MCMC
    {
        public MCMC_Gibbs():base()
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
            double cur_R0=MCMC_r0;
            
            
            //don't update, _deltaT and _duration

            /*double curLamda = lamda;
            double curQ = q;
            double curLamdaEffective = lamdaEffective;
            double curQEffective = qEffective;
            double curV0 = v0;
            */
            double sd_ka = 0.01;
            double sd_kd = 0.01;//keep it constant.
            double sd_kM = 0.00;
            double sd_conc = 0.00;
            double sd_Rmax = 0.01;//keep this unchange
            double sd_sigma = 0.01;//keep this unchange
            double sd_R0=0.01;

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
                MC_nid.setParameters(cur_ka, cur_kd, cur_kM, cur_conc, cur_Rmax, cur_R0 );
                //MC_nid.run_Attach(); 
                MC_nid.run_Detach();
                //sim_ru = MC_nid.RU_Attach;
                sim_ru = MC_nid.RU_Detach;
                //cur_loglld = logLikelihood(MC_ru_attach, sim_ru, cur_sigma);   
                cur_loglld = logLikelihood(MC_ru_detach, sim_ru, cur_sigma);
            }

            double next_conc = cur_conc;
            double next_ka = cur_ka;
            double next_kd = cur_kd;
            double next_kM = cur_kM;

            double next_Rmax = cur_Rmax ;
            double next_R0 = cur_R0;
            
            double next_sigma = cur_sigma ;
            double next_loglld=0;
            int count=1;
            while (count <= 6)
            {
                switch(count)
                {
                    case 1:
                //update the parameters
                        next_conc = Math.Exp(Math.Log(cur_conc) + sd_conc * zRand.GetRandomValue(rng));
                        break;
                    case 2:
                        next_ka = Math.Exp(Math.Log(cur_ka) + sd_ka * zRand.GetRandomValue(rng));
                        break;
                    case 3:
                        next_kd = Math.Exp(Math.Log(cur_kd) + sd_kd * zRand.GetRandomValue(rng));
                        break;
                    case 4: 
                        next_kM = Math.Exp(Math.Log(cur_kM) + sd_kM * zRand.GetRandomValue(rng));
                        break;
                    case 5:

                        next_Rmax = Math.Exp(Math.Log(cur_Rmax) + sd_Rmax * zRand.GetRandomValue(rng));
                        break;
                //double next_C0 = Math.Exp(Math.Log(cur_C0) + sd_C0 * zRand.GetRandomValue(rng));
                //double next_D0 = Math.Exp(Math.Log(cur_D0) + sd_D0 * zRand.GetRandomValue(rng));

                    case 6:
                        next_sigma = Math.Exp(Math.Log(cur_sigma) + sd_sigma * zRand.GetRandomValue(rng));
                        break;
                    case 7:
                        next_R0 = Math.Exp(Math.Log(cur_R0) + sd_R0 * zRand.GetRandomValue(rng));
                        break;
                }
                //double nextSigmaDead = Math.Exp(Math.Log(curSigmaDead) + sdSDead * zRand.GetRandomValue(rng));
                MC_nid.setParameters(next_ka, next_kd, next_kM, next_conc, next_Rmax, next_R0 );
                //MC_nid.run_Attach();
                MC_nid.run_Detach();
                //sim_ru = MC_nid.RU_Attach;
                sim_ru = MC_nid.RU_Detach;
                //next_loglld = logLikelihood(MC_ru_attach, sim_ru, next_sigma);
                next_loglld = logLikelihood(MC_ru_detach, sim_ru, next_sigma);
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
                    switch(count)
                    {
                        case 1:
                            _concArr.Add(next_conc);_conc = next_conc;
                            break;
                        case 2:
                            _kaArr.Add(next_ka);_ka = next_ka;
                            break;
                        case 3:
                            _kdArr.Add(next_kd);_kd = next_kd;
                            break;
                        case 4:
                            _kMArr.Add(next_kM);_kM = next_kM;
                            break;
                        case 5:
                            _RmaxArr.Add(next_Rmax);_Rmax = next_Rmax;
                            break;
                        case 6:
                            _sigmaArr.Add(next_sigma);_sigma = next_sigma;
                            break;
                        case 7:
                            this.MCMC_r0Arr .Add(next_R0); MCMC_r0 = next_R0;
                            break;                       
                    }
                    _cur_LogLld = next_loglld;

                }
                else
                {
                    switch (count)
                    {
                        case 1:
                            _concArr.Add(cur_conc);
                            break;
                        case 2:
                            _kaArr.Add(cur_ka);
                            break;
                        case 3:
                            _kdArr.Add(cur_kd);
                            break;
                        case 4:
                            _kMArr.Add(cur_kM);
                            break;
                        case 5:
                            _RmaxArr.Add(cur_Rmax);
                            break;
                        case 6:
                            _sigmaArr.Add(cur_sigma);
                            break;
                        case 7:
                            MCMC_r0Arr.Add(cur_R0);
                            break;    
                    }
                    _cur_LogLld = cur_loglld;
                }
                _curLldArr.Add(cur_loglld);
                _nextLldArr.Add(next_loglld);
                count++;

                
            }//end of while
            //save next paramters too
            writerNextParameter.WriteLine(steps + "\t" + next_ka + "\t" + next_kd + "\t" + next_kM + "\t"
                        + next_conc + "\t" + next_Rmax + "\t" + next_sigma + "\t" + next_R0 + "\t" + cur_loglld + "\t" + next_loglld);

            if (steps == MC_total_cycles - 1)
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
