using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using Meta.Numerics.Statistics.Distributions;
using System.IO;


namespace BayesianEstimateLib
{
    /// <summary>
    /// the class doing the Bayesian Estimation.
    /// </summary>
    public abstract class MCMC
    {
        /// <summary>
        /// empty contructor
        /// </summary>
        public MCMC()
        {
            //initialize the arrays
            _kaArr = new List<double>();
            _kdArr = new List<double>();
            _kMArr = new List<double>();
            _concArr = new List<double>();
            _RmaxArr = new List<double>();
            _sigmaArr = new List<double>();
            _curLldArr = new List<double>();
            _nextLldArr = new List<double>();
            MCMC_r0Arr = new List<double>();
            

            MC_time_attach = new List<double>();
            MC_ru_attach = new List<double>();

            MC_time_detach=new List<double>();
            MC_ru_detach = new List<double>();


            zRand = new NormalDistribution();
            rng = new Random(AccessoryLib.AceessoryLib.SEED);

            uRand = new UniformDistribution();
            rng2 = new Random(AccessoryLib.AceessoryLib.SEED+1);

            uRand2 = new UniformDistribution();
            rng3 = new Random(AccessoryLib.AceessoryLib.SEED+2);

            //initialize the values
            _deltaT = 0.01;
            //
            nextParameters = new Dictionary<string,double>();
        }

        public void getInput(List<double> _tArr_attach, List<double> _rArr_attach, List<double> _tArr_detach, List<double> _rArr_detach)
        {
            MC_time_attach = _tArr_attach;
            MC_ru_attach = _rArr_attach;
            _duration_attach = MC_time_attach.Max()+0.0;

            MC_time_detach = _tArr_detach;
            MC_ru_detach = _rArr_detach;
            _duration_detach = MC_time_detach.Max() + 0.0;
        }
       
        public void run()
        {
            //first run simulation to get cycle data
            //set up the initial values
            _ka = _ka0;
            _kd = _kd0;
            _kM = _kM0;
            _conc = _conc0;
            _Rmax = _Rmax0;
            _sigma = _sigma0;
            MCMC_r0 = MCMC_r00;

            //MC_nid = new Langmuir( _ka, _kd, _conc, _Rmax, MCMC_r0, _duration_attach, _duration_detach, _deltaT);

            MC_nid = new NumericalIntegrationOfDynamics(_ka, _kd, _conc, _Rmax,MCMC_r0, _duration_attach,_duration_detach, _kM, _deltaT );

            //MC_nid = new MassTransportModel (_ka, _kd, _conc, _Rmax,MCMC_r0, _duration_attach,_duration_detach, _kM, _deltaT );
            Console.WriteLine("start updating the input time arrays.............");
            
            
            
            MC_nid.addValueSToTimeArr_attach(MC_time_attach);
            Console.WriteLine("done for attaching...........");

            MC_nid.addValueSToTimeArr_detach(MC_time_detach);
            Console.WriteLine("done for detaching...........");
            
//MC_total_cycles =50000;
            _cur_LogLld = 1E-200;
            //now go through 
            Console.WriteLine(">>>>>Starting.........\r\n");
            Console.Out.Flush();
            this.writerNextParameter = new StreamWriter("simulation_Nextset.txt");
            writerNextParameter.WriteLine("line\tka\tkd\tkM\tconc\tRmax\tSigma\tcurLLD\tnextLLD");
            for (int i = 0; i < MC_total_cycles; i++)
            {
                MCMCStep(_cur_LogLld,i);
                if (i % 500 == 0)
                {
                    Console.WriteLine("loop:" + i + "\r\n");
                    Console.Out.Flush();
                }
            }
            writerNextParameter.Close();
            Console.WriteLine("Done.............\r\n");
            Console.Out.Flush();
        }
        /// <summary>
        /// to do one step of MCMC MH algorithm
        /// </summary>
        protected abstract void MCMCStep(double prevLogLLD, int steps);
       

        public void writeOutput()
        {
            StreamWriter writer = new StreamWriter("MCMC_run.txt");
            writer.WriteLine("line\tka\tkd\tkM\tconc\tRmax\tSigma\tR0\tcurLLD\tnextLLD");
            int count = 0;
            for (int i = 0; i < _kaArr.Count; i++)
            {
                writer.WriteLine(i+"\t"+_kaArr[i] + "\t" + _kdArr[i]+ "\t" +_kMArr[i]+"\t"
                    + _concArr[i] + "\t" + _RmaxArr[i] + "\t" + _sigmaArr[i] + "\t" + MCMC_r0Arr[i] + "\t" + _curLldArr[i] + "\t" + _nextLldArr[i]);
                if (i > MC_burn_in)
                {
                    count++;
                    mean_conc += _concArr[i];
                    mean_ka += _kaArr[i];
                    mean_kd += _kdArr[i];
                    mean_kM += _kMArr[i];
                    mean_Rmax += _RmaxArr[i];
                    mean_sigma += _sigmaArr[i];
                    mean_r0 += this.MCMC_r0Arr[i];
                }
            }

            mean_conc /= count;// (MC_total_cycles - MC_burn_in);
            mean_ka /= count;//(MC_total_cycles - MC_burn_in);
            mean_kd /= count;//(MC_total_cycles - MC_burn_in);
            mean_kM /= count;//(MC_total_cycles - MC_burn_in);
            mean_Rmax /= count;//(MC_total_cycles - MC_burn_in);
            mean_sigma /= count;//(MC_total_cycles - MC_burn_in);
            mean_r0 /= count;
            writer.Close();

            Console.WriteLine("mean_conc<-" + mean_conc);
            Console.WriteLine("mean_ka<-" + mean_ka);
            Console.WriteLine("mean_kd<-" + mean_kd);
            Console.WriteLine("mean_kM<-" + mean_kM);
            Console.WriteLine("mean_sigma<-" + mean_sigma);
            Console.WriteLine("mean_r0<-" + mean_r0);
            Console.WriteLine("mean_Rmax<-" + mean_Rmax);
        }
        protected double logLikelihood(List<double> _obs, List<double> _exp, double _sigma)
        {
            double logLL = 0; 
            //logLL += weight*logLikelihoodEventTime(firstGenDivTime, firstGenDeathTime, subSequentGenDivTime, subSequentGenDeathTime);
            for(int i=0;i<_obs.Count;i++)
            {
                logLL += logLikelihoodEach(_obs[i], _exp[i], _sigma);
            }
            return logLL;
        }

        protected double logLikelihoodEach(double _obs_i, double _exp_i, double _sigma)
        {
            double logLL;
            logLL = -0.5 * (Math.Log(2 * Math.PI * _sigma * _sigma) +
                    (_obs_i - _exp_i) * (_obs_i - _exp_i) / (_sigma * _sigma));

            return logLL;
        }

        //declare variables.
        public double _ka;
        public double _kd;
        public double _conc;
        public double _Rmax;
        public double _kM;
        public double _sigma;
        protected double _cur_LogLld;
        
//protected double _C0_detach=1;
        //protected double _D0_detach=1;
        public double MCMC_r0;//for detaching starting point 

        //starting values
        protected double _ka0=2E6;
        protected double _kd0 = 1E-4;
        protected double _conc0 = 3E-7;
        protected double _Rmax0 = 60.5;
        protected double _kM0 = 3.15E7;
        protected double _sigma0 = 1.5;
        protected double MCMC_r00=31.5;

        protected int MC_total_cycles = 30000;
        protected int MC_burn_in = 10000;

        public int MC_TotalCycles
        {
            get { return MC_total_cycles; }
            set { MC_total_cycles = value; }
        }
        public double R0
        {
            get { return MCMC_r0; }
            set { MCMC_r0 = value; }
        }

        public double R00
        {
            get { return MCMC_r00; }
            set { MCMC_r00 = value; }
        }

        public double mean_ka=0;
        public double mean_kd=0;
        public double mean_conc=0;
        public double mean_Rmax=0;
        public double mean_kM=0;
        public double mean_sigma=0;
        public double mean_r0 = 0;
        //list for remembering sequence
        protected List<double> _kaArr;
        protected List<double> _kdArr;
        protected List<double> _concArr;
        protected List<double> _RmaxArr;
        protected List<double> MCMC_r0Arr;
        protected double _duration_attach;
        protected double _duration_detach;
        protected List<double> _kMArr;
        protected double _deltaT;
        protected List<double> _sigmaArr;
        protected List<double> _curLldArr;
        protected List<double> _nextLldArr;

        //data 
        protected List<double> MC_time_attach;
        protected List<double> MC_ru_attach;

        protected List<double> MC_time_detach;
        protected List<double> MC_ru_detach;

        public SimulationSPR  MC_nid;

        //random generator
        protected NormalDistribution zRand;
        protected Random rng;

        protected UniformDistribution uRand;
        protected Random rng2;

        protected UniformDistribution uRand2;
        protected Random rng3;

        protected StreamWriter writerNextParameter;

        public Dictionary<string, double> nextParameters;
    }//end of class
}
