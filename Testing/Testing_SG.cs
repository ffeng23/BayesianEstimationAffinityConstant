using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using Savitzky_GolaySmoothingLib;
using BayesianEstimateLib;

using AccessoryLib;
namespace Testing
{
    class Testing_SG
    {
        static void Main(string[] args)
        {
            //first testing factorization simplying
            List<BigInteger> temp=new List<BigInteger>{-4,2};
            BigInteger num_r = -4*34*17*24;
            BigInteger den_r = 2 * 34 * 17;
            Console.WriteLine("the input number:" + num_r + "/" + den_r);
            SG_Smoothing.SimplifyNumDenom(ref num_r, ref den_r);

            Console.WriteLine("the returning list:" + num_r + "/" + den_r);

            Console.WriteLine("********************************");
            Console.WriteLine("doing general factorial...............");
            BigInteger tempInt = SG_Smoothing.GeneralFactorial(2*20+5+1, 5+1);
            //Console.WriteLine("factorial of 10 and 9 is " + tempInt);

            Console.WriteLine("********************************");
            Console.WriteLine("doing grand polynomial...............");
            int index=1;int orderOfDerivative=1;int degreeOfPolynomial=1;int frameSize=3;
            temp=SG_Smoothing.GrandPolynomial( index, orderOfDerivative, degreeOfPolynomial, frameSize);

            Console.WriteLine("the grand polynomial is i="+index+";orderOfDerivative="+orderOfDerivative+
                ";degreeOfPolynomial="+degreeOfPolynomial+";frameSize="+frameSize + "is "+temp[0]+"/"+temp[1]);

            Console.WriteLine("********************************");
            Console.WriteLine("doing weight...............");
            int positionOfpoints = -10;//this is the t in the reference paper
            index = 0; frameSize =10; degreeOfPolynomial = 1; orderOfDerivative = 0;
            List<BigInteger> tempInt64 = SG_Smoothing.Weight(index, positionOfpoints, frameSize, degreeOfPolynomial, orderOfDerivative);
            Console.WriteLine("the weight coefficient is i=" + index + ";orderOfDerivative=" + orderOfDerivative +
                ";degreeOfPolynomial=" + degreeOfPolynomial + ";frameSize=" + frameSize + "is " + tempInt64[0] + "/" + tempInt64[1]);
            SG_Smoothing.SimplifyNumDenom(temp);
            Console.WriteLine("the weight is i=" + index + ";orderOfDerivative=" + orderOfDerivative +
                ";degreeOfPolynomial=" + degreeOfPolynomial + ";frameSize=" + frameSize + "is " + tempInt64[0] + "/" + tempInt64[1]);

            Console.WriteLine("***************testing matrix*****************");
            Console.WriteLine("doing matrix...............");
            frameSize=20;degreeOfPolynomial=2;
            SG_Smoothing sgs = new SG_Smoothing(frameSize,degreeOfPolynomial );
            /*Dictionary<int, List<List<BigInteger>>> cfm=sgs.GenerateCoefficients(0);
            Console.WriteLine("total number of matrix:" + cfm.Count);
            for (int t = -1*frameSize ; t <= frameSize; t++)
            {
                List<List<BigInteger>> cfm_t = cfm[t];
                Console.Write("t=" + t + ":");
                for (int i = 0; i < cfm_t.Count; i++)
                {
                    Console.Write(cfm_t[i][0] +"/"+cfm_t[i][1]+ "\t");
                }
                Console.Write("\n");
            }

            */
            Console.WriteLine("***************testing smothing*****************");
            Console.WriteLine("doing IO...............");
            Dictionary<int, List<Double>> dt=DataIO.ReadDataTable("simulation_detach_noise.txt");
            Console.WriteLine("there are " + dt.Count + " columns and " + dt[0].Count + " rows.");
            Dictionary<int, List<Double>> dt_short = new Dictionary<int, List<double>>();
            List<int> key=dt.Keys.ToList();
            dt_short.Add(dt.Keys.ToList()[0], new List<double>());
            dt_short.Add(dt.Keys.ToList()[1], new List<double>());
            for (int i = 0; i < dt[0].Count; i=i+900)
            {
                dt_short[0].Add(dt[0][i]);
                dt_short[1].Add(dt[1][i]);
            }

            //now have a input array for doing smoothing
            DataIO.WriteDataTable("simulation_detach_noise_short.txt", dt_short);

            //nowing doing the smoothing
            List<double> dt_smoothed = sgs.Smooth(0, dt_short[1]);
            List<string> header=new List<string>{"time","RU"};
            DataIO.WriteDataTable(dt_short[0], dt_smoothed, "simulation_detach_noise_smoothed.txt", header);

            //frameSize = 50; degreeOfPolynomial = 1;
            //sgs.SetParameters(frameSize, degreeOfPolynomial);
            //List < double> dt_smoothed2 = sgs.Smooth(0, dt_smoothed);
            //DataIO.WriteDataTable(dt_short[0], dt_smoothed2, "simulation_detach_noise_smoothed2.txt", header);
            //nowing doing the slope
            List<double> dt_slope = sgs.Smooth(1, dt_short[1]);
            List<string> header_slope = new List<string> { "time", "RU" };
            DataIO.WriteDataTable(dt_short[0], dt_slope, "simulation_detach_noise_slope.txt", header_slope);
        }//main function
    }//end of class
}//end of namespace
