using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NelderMeadMethod;
using BayesianEstimateLib;
using System.IO;

using AccessoryLib;
using System.Numerics;
using Models;

namespace testingNelderMead1D
{
    //public delegate double LogDistributionFuctionDelegate(double x, double functionNormConstant);
    class testingNelderMead1D
    {
        static void Main(string[] args)
        {
            Console.WriteLine("start testing logNormal distribution..........");
            //first testing the logNormalDist
            int count=1000;
            List<double> x = new List<double>(count), y=new List<double>(count);
            double low=-250;
            for (int i = 0; i < count; i++)
            {
                x.Add(low + i * 0.5);
                y.Add(logNormalDist(x[i], 0));
            }

            //now write the output
            Console.WriteLine("Write the output............");
            DataIO.WriteDataTable(x, y, "normalDist.txt",  new List<string> { "x", "y" });
            Console.WriteLine("Done............");

            double funcNormalConst = 0;
            Console.WriteLine("Start testing nelder mead algorithm..........");
            double nonzero = NelderMead1D.FindNonZeroValue(logNormalDist, 0, 3000, 0, 3000, NelderMeadMethod.NelderMead1D.LOG_LIMIT, ref funcNormalConst);
            Console.WriteLine("the nonzero value is " + nonzero);
            Console.WriteLine("doing............");

            //testing the conditionals
            Console.WriteLine("start testing nelder mead algorithm..........");

            //testing the solve cubic equation
            bool realRoot = false;
            List<Complex> roots=AccessoryLib.AccessoryLib.SolveCubic(1,-6, 11, -6, out realRoot );
            if (realRoot)
            {
                Console.WriteLine("there are real number roots");
            }
            foreach (Complex c in roots)
            {
                Console.WriteLine(c.ToString());
            }

            //*******************************now start doing testing of Gibbs sampler on two state model of spr
            Console.WriteLine("Hello world!!!************doing the gibbs on two state model");

            //Start doing the 
            FitController_UnifiedTwoState fc = new FitController_UnifiedTwoState(1E-7, 1000);

            fc.SetupModel();
            List<List<double>> dist = fc.Run(1500);
            /*
            //write ing output
            Console.WriteLine("writing the output file.......");
            StreamWriter writer = new StreamWriter("learReg.txt");
            writer.WriteLine("kon_if\tkoff_if\tka_if\tkd_if\tkon_cs\tkoff_cs\tka_cs\tkd_cs\tRmax\tau");
            for (int i = 0; i < dist[0].Count; i++)
            {
                for (int j = 0; j < dist.Count; j++)
                {
                    writer.Write(dist[j][i]);
                    if (j != dist.Count - 1)
                    {
                        writer.Write("\t");
                    }
                    else
                        writer.Write("\r\n");
                }
            }
            Console.WriteLine("done!!!!!!!");
            writer.Close();*/
            //end of code.
            Console.WriteLine("done, hit enter to quit");
            Console.ReadLine();

        }


        //define a normal distribution for testing
        /// <summary>
        /// this is method for testing
        /// </summary>
        /// <param name="_x">the input for doing logNormal</param>
        /// <param name="normalConst">this is the constant to normalize the result. we actually do not use this for the actual function
        ///             We use this only for compatibility issue. since NelderMead1D.FindNonZerValue needs this form of function to work
        ///             For this function logNormal, it is easy and no need for this. But for other complicated functions, we need this
        ///             constant to normalize the output and find the value.
        /// </param>
        /// <returns></returns>
        static double logNormalDist(double _x, double normalConst)
        {
            double y;
            double mu = 2500, sigma = 0.1;

            y = Math.Log(1 / Math.Sqrt(2 * Math.PI * sigma * sigma)) - (0.5 * (_x - mu) * (_x - mu) / (sigma*sigma));
            return y;
        }

    }
}
