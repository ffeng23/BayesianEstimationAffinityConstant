using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using AdaptiveRejectionSampling;
using Meta.Numerics.Statistics.Distributions;
using NelderMeadMethod;
using System.IO;
using AccessoryLib;

namespace GibbsSampler
{
    /// <summary>
    /// this is the function used to udpate the distribution function based on the newly drawn paramters
    /// </summary>
    /// <param name="parameters">contains the paramter that is updated so far</param>
    /// <param name="index">the index for which the distribution function should be update and returns</param>
    /// <returns></returns>
    public delegate LogDistributionFuctionDelegate UpdateDistributionDelegate(List<double> parameters, int index); 
    /// <summary>
    /// finally this is the implementation of Gibbs sampler using the ARMS algorithm to draw samples from arbitrary target distribution
    /// </summary>
    public class GibbsSampler
    {
        /// <summary>
        /// empty constructor, disabled
        /// </summary>
        protected GibbsSampler()
        {
            //empty contructor, disabled
        }

        /// <summary>
        /// constructor
        /// </summary>
        /// <param name="_initials">the initial values for the parameters, this is necessary, since we are doing ARMS</param>
        ///  <param name="_updateDistribution">the function used to update the distribution function based on new parameter values</param>
        /// <param name="_bounds">the lower and upper bound for each parameter</param>
        public GibbsSampler( List<double> _initials, UpdateDistributionDelegate _updateDistribution, List<List<double>> _bounds)
        {
            this.CP_Initials = _initials;
            this.C_UpdateDistributionDelegate = _updateDistribution ;
            this.C_Bounds = _bounds;
            rngs = new List<Random>(_initials.Count);
            for(int i=0;i<_initials.Count;i++)
            {
                rngs.Add(new Random(AccessoryLib.AceessoryLib.SEED+i));
            }
        }

        public List<List<double>> Run(int _numberOfSamples)
        {
            Console.WriteLine("starting.......");
            //initalize the arrays.
            List<List<double>> samples = new List<List<double>>(this.CP_Initials.Count);
            AdaptiveRejectionMetropolisSampling arms;// = new AdaptiveRejectionMetropolisSampling>(this.CP_Initials.Count);
            //List<LogDistributionFuctionDelegate> distributions=new List<LogDistributionFuctionDelegate>(this.CP_Initials.Count);
            //update/initalize all the distribution functions
            Console.WriteLine("initializing distribution............");
            /*for(int i=0;i<this.CP_Initials.Count;i++)
            {
                distributions.Add(C_UpdateDistributionDelegate(this.CP_Initials,i));
            }*/
            Console.WriteLine("initialize bayesisna samplings.........");
            List<double> current=this.CP_Initials;
            List<double> previous = new List<double>();//this one is used to remember the previous sample in order to be used for 1d search of nonzero for nelder mead method
            for (int i = 0; i < this.CP_Initials.Count; i++)
            {
                samples.Add(new List<double>(_numberOfSamples));
                samples[i].Add(this.CP_Initials[i]);
                double lower = this.C_Bounds[i][0] ;
                if (Double.IsInfinity(lower))
                {
                    lower = -1E20;
                }
                double upper = this.C_Bounds[i][1];
                if (Double.IsInfinity(upper))
                {
                    upper = 1E20;
                }
                
                if(this.CP_Initials[i]!=(lower +upper)/2)
                    previous.Add((lower + upper) / 2);
                else
                    previous.Add((lower+this.CP_Initials[i])/2);

                //initialize the ARMS for sampling
                //arms.Add( new AdaptiveRejectionMetropolisSampling(this.CP_Initials[i], distributions[i], this.C_Bounds[i][0], this.C_Bounds[i][1]);
            }

            Console.WriteLine("Start drawing samples..............");
            //run to generatate samples
            StreamWriter writer = new StreamWriter("learReg.txt");
            //writer.WriteLine("ka\tkb\tkM\tconc\tRmax\tR0\tVar");
            int increPercent = 5;
            int runningPercent = 0;
            Console.Write("Progress:");
            for (int i = 1; i < _numberOfSamples; i++)
            {
                if (((double)i) / _numberOfSamples * 100 >= runningPercent)
                {
                    Console.Write("..." + runningPercent + "%(" + i + "/"+_numberOfSamples+")" );
                    runningPercent += increPercent;
                }
                //for each different paramter, gerenate one sample from its target distribution
                for (int j = 0; j < samples.Count; j++)
                {
                    //Console.WriteLine("\tsub " + j + " steps; seed is "+AccessoryLib.AceessoryLib.SEED );
                    LogDistributionFuctionDelegate func = C_UpdateDistributionDelegate(current, j);
                    arms = new AdaptiveRejectionMetropolisSampling( previous[j],current[j], func, this.C_Bounds[j][0], this.C_Bounds[j][1]);
                    //Console.WriteLine("\t\tfinished with ARMS set up...");
                    double newSample = arms.GetRandomSample(rngs[j]);
                    //Console.WriteLine("\t\tfinished with sample drawing..."+newSample );
                    samples[j].Add(newSample);
                    if (newSample != current[j])
                    {
                        previous[j] = current[j];
                    }//otherwise we won't update keep the previous one in case we have same values to crash the code.
                    current[j] = newSample;
                    
                    writer.Write(current[j]);
                    if (j != samples.Count - 1)
                    {
                        writer.Write("\t");
                    }
                    //now update the distribution function for the next parameter
                    //distributions[(j + 1) % (samples.Count )] = C_UpdateDistributionDelegate(current, (j + 1) % (samples.Count ));
                    //Console.WriteLine("\tfinished with sub steps...");
                    
                }//for different parameters
                writer.WriteLine(""); writer.Flush();
                

            }//for different sample rounds
            Console.WriteLine("done..........");
            writer.Close();
            return samples;
        }
        /// <summary>
        /// using MH like algorithm to estimate the initial values of parameter, the goal is to have the initials to make all the functions to be nonzero
        /// </summary>
        /// <returns></returns>
        /*
        public List<double> GenerateInitials(List<LogDistributionFuctionDelegate> _logDist)// List<DistributionFuctionDelegate> _dist)
        {
            Random rng = new Random();
            NormalDistribution zRand = new NormalDistribution();
            LogDistributionFuctionDelegate logDistribution;
            
            List<Double> currentLogDistFuncValue=new List<double>(this.CP_Initials.Count);
            for(int i=0;i<this.CP_Initials.Count;i++)
            {
                currentLogDistFuncValue.Add(-1);
            }

            bool flagAllNonZero=true;
            //the outer loop to 

            while (true)
            {
                for (int i = 0; i < this.CP_Initials.Count; i++)
                {
                    logDistribution = this.C_UpdateDistributionDelegate(this.CP_Initials, i);
                    if (Math.Exp(logDistribution(this.CP_Initials[i])) <= 0)//not good for this one
                    {
                        flagAllNonZero = false;

                        //we need to check for the where we are accepting this one or not
                    }
                    
                }
            }
        }
        */
        //**********declaration of members
        List<double> CP_Initials;
        UpdateDistributionDelegate C_UpdateDistributionDelegate;
        List<List<double>> C_Bounds;

        private List<Random> rngs;
    }//end of class
}//end of namespace
