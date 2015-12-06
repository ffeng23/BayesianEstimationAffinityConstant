using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using GibbsSampler;
using Meta.Numerics.Statistics.Distributions;
using AdaptiveRejectionSampling;
using System.IO;
//using GibbsSampler;
using BayesianEstimateLib;
using NelderMeadMethod;

namespace TestingLinearRegression
{
    /// <summary>
    /// testing linear regression using Gibbs sampler
    /// </summary>
    class TestingLinearRegression
    {
        static void Main(string[] args)
        {

            LogDistributionFuctionDelegate fd = new LogDistributionFuctionDelegate(dist2);
            Console.WriteLine("First try:" + fd(5,0));
                fd=UpdateDistribution(new List<double>() { 1, 2 }, 1);
            Console.WriteLine("Second try:" + fd(5,0));

            Console.WriteLine("Starting doing testing linear regression...................");
            LinearRegression lr = new LinearRegression(3000, 1500, 1.5);
            int numberSamples = 40;
            List<double> x = new List<double>(numberSamples);

            for (int i = 0; i < numberSamples; i++)
            {
                x.Add(i);
            }

            lr.GetYValues(x);
            List<double> y = lr.C_Y;

            StreamWriter writer = new StreamWriter("sampleLinear.txt");
            writer.WriteLine("X\tY");
            
            for (int i = 0; i < numberSamples; i++)
            {
                
                writer.WriteLine(x[i]+ "\t" + y[i]);
                //Console.Write("," + x_array[i]);
            }
            writer.Close();

            //now need to fit the model
            Console.WriteLine("*****testing. fitting the model......");
            List<List<double>> bounds = new List<List<double>>();
            bounds.Add(new List<double> { Double.NegativeInfinity , Double.PositiveInfinity  });
            bounds.Add(new List<double> { Double.NegativeInfinity, Double.PositiveInfinity });
            bounds.Add(new List<double> { 0,100 });
            Console.WriteLine("building the gibbs sampler...........");
            List<int> lstFunc=new List<int>();
            lstFunc.Add(0 /*slope*/);
            lstFunc.Add(1 /*intercept*/);
            lstFunc.Add(2/*var*/);
            lr.setFunctionDelegateForUpdating(lstFunc);
            GibbsSampler.GibbsSampler gb = new GibbsSampler.GibbsSampler(new List<double> { 6000, 100, 10 /* */}, lr.updateFunctionDistribution, bounds);
            Console.WriteLine("Done...........");
            
            List<List<double>> output = new List<List<double>>();

            //testing search for nozero new method
            //PiecewiseUniformWithExponentialTails pw = new PiecewiseUniformWithExponentialTails(1000, Double.NegativeInfinity, Double.PositiveInfinity, lr.SlopeCondistionalDist);
            //lr.updateFunctionDistribution(new List<double> { 1, 5, 1 }, 1);
            Console.WriteLine("Seraching using new............");
            //double nozero=pw.SearchNonZeroPointNew(lr.SlopeCondistionalDist);

            //Console.WriteLine("non zero value is " + nozero);


            output=gb.Run(1500);

            //write ing output
            Console.WriteLine("writing the output file.......");
            writer = new StreamWriter("learReg.txt");
            writer.WriteLine("a\tb\tau");
            for(int i=0;i<output[0].Count ;i++)
            {
                for(int j=0;j<output.Count;j++)
                {
                    writer.Write(output[j][i]);
                    if (j != output.Count-1)
                    {
                        writer.Write("\t");
                    }
                    else
                        writer.Write("\r\n");
                }
            }
            Console.WriteLine("done!!!!!!!");
            writer.Close();


        }//end of main

        //function to update the distribution function based on the new parameter
        static LogDistributionFuctionDelegate UpdateDistribution(List<double> _params, int index)
        {
            a = _params[0];
            b = _params[1];
            if (index == 0)
                return dist1;
            else
                return dist2;
        }

        static double dist1(double _x, double _x2)
        {
            //a = 0;

            return _x*a + 1;
        }
        static double dist2(double _x, double _x2)
        {
            //b= 0;

            return _x / b + 1;
        }
        static double a=0;
        static double b=1;
    }
}
