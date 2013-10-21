using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

using AdaptiveRejectionSampling;
using NelderMeadMethod;

namespace Testing2
{
    class Program
    {
        static void Main(string[] args)
        {
            double x = Double.PositiveInfinity;
            Console.WriteLine("log(Inf):"+normalDist(x,0));
            Console.WriteLine("log(0):" + normalDist(0,0));

            PiecewiseUniformWithExponentialTails puwt = new PiecewiseUniformWithExponentialTails(11,  -10  , Double.PositiveInfinity, new LogDistributionFuctionDelegate(normalDist),0,100);

            //double starting=puwt.SearchNonZeroPoint(normalDist, -1E6,Double.PositiveInfinity);
            //Console.WriteLine("the starting is " + starting);
            

            List<double> supportpoints = puwt.SupportPoints;

            Console.WriteLine("total number1 of points is " + supportpoints.Count);
            Console.Write("s={");
            for (int i = 0; i < supportpoints.Count; i++)
            {
                Console.Write( supportpoints[i] + ",");
            }
            Console.WriteLine("}");

            //now testing the envelope function of log target
            List<List<double>> envLogTarget = puwt.EnvelopeFunctionOfLogTarget;
            Console.WriteLine("the total pieces in the envelope function is " + envLogTarget.Count);

            for (int i = 0; i < envLogTarget.Count; i++)
            {
                Console.WriteLine("i=" + i + ":" + envLogTarget[i][0] + "x+" + envLogTarget[i][1]);
            }

            //***now testing the pdf cdf etc
            Console.WriteLine("total pdf normalization factor is " + puwt.NormalizationFactorOfProb);

            Console.Write("piecewisePDF: {");
            for (int i = 0; i < puwt.PiecewisePDF.Count; i++)
            {
                Console.Write("," + puwt.PiecewisePDF[i]);
            }
            Console.WriteLine("}");

            Console.Write("piecewiseCDF: {");
            for (int i = 0; i < puwt.PiecewiseCDF.Count; i++)
            {
                Console.Write("," + puwt.PiecewiseCDF[i]);
            }
            Console.WriteLine("}");

            //**********testing insertion of new points
            puwt.AddOnePointToSupportPointsArray(19.232335);
            puwt.AddOnePointToSupportPointsArray(0.5);
            puwt.AddOnePointToSupportPointsArray(-0.5);
            Console.WriteLine("***********add one more point************");

            Console.WriteLine("total number1 of points is " + supportpoints.Count);
            Console.Write("s={");
            for (int i = 0; i < supportpoints.Count; i++)
            {
                Console.Write(supportpoints[i] + ",");
            }
            Console.WriteLine("}");

            envLogTarget = puwt.EnvelopeFunctionOfLogTarget;

            Console.WriteLine("the total pieces in the envelope function is " + envLogTarget.Count);

            for (int i = 0; i < envLogTarget.Count; i++)
            {
                Console.WriteLine("i=" + i + ":" + envLogTarget[i][0] + "x+" + envLogTarget[i][1]);
            }


            Console.WriteLine("total pdf normalization factor is " + puwt.NormalizationFactorOfProb);

            Console.Write("piecewisePDF: {");
            for (int i = 0; i < puwt.PiecewisePDF.Count; i++)
            {
                Console.Write("," + puwt.PiecewisePDF[i]);
            }
            Console.WriteLine("}");

            Console.Write("piecewiseCDF: {");
            for (int i = 0; i < puwt.PiecewiseCDF.Count; i++)
            {
                Console.Write("," + puwt.PiecewiseCDF[i]);
            }
            Console.WriteLine("}");

            Console.WriteLine("testing.......sampling..........of piecewise function*******************");
            Random rng = new Random(AccessoryLib.AceessoryLib.SEED );
            double temp = puwt.GetRandomValue(rng);
            Console.WriteLine("the sample is "+temp);

            Console.WriteLine("expHnOfX is " + puwt.ExpHnOfX(temp));

            Console.WriteLine("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
            Console.WriteLine("Starting testing ARMS..............");
            AdaptiveRejectionMetropolisSampling arms = new AdaptiveRejectionMetropolisSampling(0.5,0.6, new LogDistributionFuctionDelegate(betaDist),0,1);
            List<double> x_array = new List<double>();
            StreamWriter writer = new StreamWriter("sample.txt");
            writer.WriteLine("id\tsample");
            //Console.Write("{");
            for (int i = 0; i < 10000; i++)
            {
                x_array.Add(arms.GetRandomSample(rng));
                writer.WriteLine(i + 1 + "\t" + x_array[i]);
                //Console.Write("," + x_array[i]);
            }
            writer.Close();
            //Console.WriteLine("}");
        }//end of main

        static double normalDist(double _x, double _normalizationConstant)
        {
            double sigma=0.5;
            double mu=0;
            return -1*(Math.Log( Math.Sqrt(2 * Math.PI * sigma*sigma))) +(-0.5*(_x-mu)*(_x-mu) / (sigma * sigma));
            
        }

        static double betaDist(double _x, double _nomalizationConstant)
        {
            double alpha = 1;
            double beta = 3;
            double constant = 10;
            return Math.Log(1/constant*Math.Pow(_x, alpha-1)*Math.Pow((1-_x),beta-1));

        }
    }//end of class
}//end of namespace
