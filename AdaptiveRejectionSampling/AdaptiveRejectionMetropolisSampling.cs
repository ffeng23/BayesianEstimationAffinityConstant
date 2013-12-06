using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NelderMeadMethod;
using AccessoryLib;

namespace AdaptiveRejectionSampling
{
    
    /// <summary>
    /// this is class implementing Adaptive Rejection Metropolis Sampling. In this class, it doesn't assume a log-concave target function
    /// This implementation is mainly based on "Adative Rejection Metropolis Sampling within Gibbs Sampling. Gilks, WR, Best, NG and Tan, KKC
    /// 1995. Journal of the Royal Statistical Society.Series C (Applied Statistics). Vol. 44. No. 4 pp. 455
    /// </summary>
    public class AdaptiveRejectionMetropolisSampling
    {
        private AdaptiveRejectionMetropolisSampling()
        {
            //empty contructor
            //disabled for now.
        }
        /// <summary>
        /// constructor, the lower or upper bound is not Very important, since we will search for nonzero upper or lower bound anyway, so if we don't know the bounds, we
        /// can set it up as infinity. if we set it as noninfinity, they have to be reasonable, at least lead to non zero value of the distribution, because we assume this 
        /// happens in the piecewise envelope function.
        /// </summary>
        /// <param name="_dist">the target distribution we tried to sample from</param>
        /// <param name="_lowerBound">lower bound of the distribution</param>
        /// <param name="_upperBound">upper bound of the distribution</param>
        /// <param name="_numberOfSupportPoints">how many points we start with. normally this will create evenly spaced intervals</param>
        public AdaptiveRejectionMetropolisSampling(Double _X_initialPrevious, Double _X_initialCurrent, LogDistributionFuctionDelegate _LogDist, double _lowerBound=Double.NegativeInfinity , double _upperBound=Double.PositiveInfinity , int _numberOfSupportPoints=20)
        {
            //now get the proposal distribution
            this.CP_ProposalDistribution = new PiecewiseUniformWithExponentialTails(_numberOfSupportPoints, _lowerBound, _upperBound, _LogDist, _X_initialPrevious, _X_initialCurrent );
            this.CP_LogTargetDistribution = _LogDist;
            

            this.CP_uniformRng1 = new Random(AccessoryLib.AceessoryLib.SEED);
            this.CP_uniformRng2 = new Random(AccessoryLib.AceessoryLib.SEED+1);

            this.CP_X_cur = _X_initialCurrent;
        }

        public double GetRandomSample(Random rng)
        {
            //by now we should have had the support points initialized and envelope function ready, 
            //now we need to run the ARMS algorithm to draw one sample
            double u, X, fX, exphX, XA, XM;
            while(true)
            {
                u = CP_uniformRng1.NextDouble(); //step 2
                X = CP_ProposalDistribution.GetRandomValue(rng); //step 1
                fX=Math.Exp(CP_LogTargetDistribution(X, this.CP_ProposalDistribution.FunctionNormConstant ));
                
                exphX=this.CP_ProposalDistribution.ExpHnOfX(X);
                //Console.WriteLine("calling in loop fx:"+fX+"; expHx:"+exphX+";u:"+u);

                //step 3
                if (u > fX / exphX)
                {
                    //need to add this point to the support point array
                    this.CP_ProposalDistribution.AddOnePointToSupportPointsArray(X);
                    continue;
                }
                else
                {
                    XA = X;
                    break;
                }
            }

            //step 4, MH step
            u = CP_uniformRng2.NextDouble();

            //step 5,
            double f_XA=fX;//evaluation of XA is done in above, we don't have to do it again.
            double f_Xcur=Math.Exp(this.CP_LogTargetDistribution(this.CP_X_cur, this.CP_ProposalDistribution.FunctionNormConstant));
            double expHnXA=exphX;//evaluation of expHnX is done above, we don't have to do it again.
            //Console.WriteLine("calling in MH step");
            double expHnXcur=this.CP_ProposalDistribution.ExpHnOfX(this.CP_X_cur);
            
            double r1, r2;
            //take min(f(Xcur), expHn(Xcur)) 
            if(f_Xcur<expHnXcur)
            {
                r1=f_XA*f_Xcur;
            }
            else
            {
                r1=f_XA*expHnXcur;
            }
            //take min(f(XA), expHn(XA))
            if(f_XA<expHnXA)
            {
                r2=f_Xcur*f_XA;
            }
            else
            {
                r2=f_Xcur*expHnXA ;
            }

            if(u>r1/r2)//rejection
            {
                XM = CP_X_cur;
            }
            else
            {
                XM=XA;
            }

            this.CP_X_cur = XM;
            return XM;
        }


        //**************memeber declaration
        PiecewiseUniformWithExponentialTails CP_ProposalDistribution;
        LogDistributionFuctionDelegate CP_LogTargetDistribution;
        Random CP_uniformRng1;
        Random CP_uniformRng2;
        Double CP_X_cur;
    }//end of class
}
