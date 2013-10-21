using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NelderMeadMethod;
using BayesianEstimateLib;

namespace testingNelderMead1D
{
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
                y.Add(logNormalDist(x[i]));
            }

            //now write the output
            Console.WriteLine("write the output............");
            DataIO.WriteDataTable(x, y, "normalDist.txt",  new List<string> { "x", "y" });
            Console.WriteLine("Done............");

            Console.WriteLine("start testing nelder mead algorithm..........");
            double nonzero=NelderMead1D.FindNonZeroValue(logNormalDist, 0, 3000, 0, 3000, -5);
            Console.WriteLine("the nonzero value is " + nonzero);
            Console.WriteLine("doing............");

            //testing the conditionals
            Console.WriteLine("start testing nelder mead algorithm..........");

        }


        //define a normal distribution for testing

        static double logNormalDist(double _x)
        {
            double y;
            double mu = 2500, sigma = 0.1;

            y = Math.Log(1 / Math.Sqrt(2 * Math.PI * sigma * sigma)) - (0.5 * (_x - mu) * (_x - mu) / (sigma*sigma));
            return y;
        }

    }
}
