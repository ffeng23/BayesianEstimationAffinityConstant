using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BayesianEstimateLib;
using GibbsSampler;
using AdaptiveRejectionSampling;
using System.IO;
using NelderMeadMethod;
using Models;

using AccessoryLib;
namespace SPR_GibbsSampler
{
    class SPRGibbs
    {
        static void Main(string[] args)
        {
            /*_ka = 1.1E5;
            _kd = 1.28E-4;
            _conc = 4E-8;
            _Rmax = 1520;

            _duration_attach = 400;
            _duration_detach = 400;

            _kM = 3.15E8; //this is now in RU/S
            _deltaT = 0.01;
            _R0 = 97.15;
            double var = 1.5;
            */
            /*_ka = 1.38E6;
            _kd = 4.8E-2;
            _conc = 0.192E-9;
            _Rmax = 481;
            */
            _ka = 1.47E5;
            _kd = 3.9E-2;
            _conc = 10E-6;
            _Rmax =120;
            
            _duration_attach = 400;
            _duration_detach = 400;

            _kM = 8E8; //this is now in RU/S
            _deltaT = 0.01;
            _R0 = 100.1;
            double var = 1.5;
            
            Console.WriteLine("run simulation....................");

            nid = new NumericalIntegrationOfDynamics(_ka, _kd, _conc, _Rmax, _R0, _duration_attach, _duration_detach, _kM, _deltaT);
            double r0 = 0;
            List<string> header = new List<string>(2);
            header.Add("time");
            header.Add("RU");
            //doing the simulation for different Rmax
            
            nid.run_Attach();
            nid.R0 = nid.RU_Attach[nid.RU_Attach.Count - 1];
            _R0 = nid.R0;
            nid.run_Detach();


            Console.WriteLine("write output..................");
            writeOutput(nid.Time_Attach, nid.RU_Attach, "simulation_attach.txt", header);
            writeOutput(nid.Time_Detach, nid.RU_Detach, "simulation_detach.txt", header);


            //writeOutput(nid.Time_Attach, nid.RU_Attach, "simulation_attach_Rmax"+(_Rmax-stepOfRmax*i) +".txt", header);
            //writeOutput(nid.Time_Detach, nid.RU_Detach, "simulation_detach_Rmax" + (_Rmax - stepOfRmax*i) + ".txt", header);

            Console.WriteLine("add noise...............");
            r0 = nid.RU_Attach[nid.RU_Attach.Count - 1];
            nid.addNoise(Math.Sqrt(var));
            
            writeOutput(nid.Time_Attach, nid.RU_Attach, "simulation_attach_noise.txt", header);
            writeOutput(nid.Time_Detach, nid.RU_Detach, "simulation_detach_noise.txt", header);

            //nid.setParameters(_ka, _kd, _kM, _conc, _Rmax, nid.R0);
            Console.WriteLine("done............");

            //now we need to shorten the data input
            List<double> shortT_A = new List<double>();
            List<double> shortR_A = new List<double>();
            List<double> shortT_D = new List<double>();
            List<double> shortR_D = new List<double>();
            for (int i = 0; i < nid.Time_Attach.Count; i=i+30)
            {
                shortT_A.Add(nid.Time_Attach[i]);
                shortR_A.Add(nid.RU_Attach[i]);
            }
            for (int i = 0; i < nid.Time_Detach.Count; i = i + 30)
            {
                shortT_D.Add(nid.Time_Detach[i]);
                shortR_D.Add(nid.RU_Detach[i]);
            }
            if (args.Count() != 2)
            {
                Console.WriteLine("writing the new output!!!!!!!!!!!!!!!!!!!!");
                writeOutput(shortT_A, shortR_A, "simulation_attach_noiseShort.txt", header);
                writeOutput(shortT_D, shortR_D, "simulation_detach_noiseShort.txt", header);
            }
            //now we have the input and testing add new values to time array and also testing the extract RUs based on the time Array
            /*Console.WriteLine("testing add new values to the time array");
            nid.addValueSToTimeArr_attach(new List<double> { 0.03, 1.03, 1.05, 2.06, 2.07 });

            nid.run_Attach();
            writeOutput(nid.Time_Attach, nid.RU_Attach, "simulation_attach.txt", header);

            nid.addValueSToTimeArr_detach(new List<double> { 0.03, 1.03, 1.05, 2.06, 2.07 });

            nid.run_Detach();
            writeOutput(nid.Time_Detach, nid.RU_Detach, "simulation_detach.txt", header);

            //now get the information
            List<double> ruExt= nid.GetRUValuesAttach(new List<double> { 0.03, 1.03, 1.05, 2.06, 2.07 });
            writeOutput(new List<double> { 0.03, 1.03, 1.05, 2.06, 2.07 }, ruExt, "simulation_attachExt.txt", header);
            ruExt = nid.GetRUValuesDetach(new List<double> { 0.03, 1.03, 1.05, 2.06, 2.07 });
            writeOutput(new List<double> { 0.03, 1.03, 1.05, 2.06, 2.07 }, ruExt, "simulation_detachExt.txt", header);
            Console.WriteLine("Done.........!");*/

            //testing read in data
            Dictionary<string, List<double> > inputDataShort;//=BayesianEstimateLib.DataIO.ReadDataTable("simulation_attach_noiseShort.txt", true, '\t', 0);


            //***********************************************GIBBS SAMPLER********************
            //now testing Gibbs sampler
            Console.WriteLine("*****testing. fitting the model......");

            //look to see whether we want to use the new simulated input or the one from the file
            string fileAttach = "simulation_attach_noiseShort.txt";
            string fileDetach = "simulation_detach_noiseShort.txt";
            if (args.Count() == 2)
            {
                Console.WriteLine("*****reading the input from disk*****");
                Console.WriteLine("reading \"" + args[0] + "\"........");
                fileAttach=args[0];
            }
                //we need to read in the input
                inputDataShort = DataIO.ReadDataTable(fileAttach, true, '\t', 0);
                shortT_A = inputDataShort["time"];
                shortR_A = inputDataShort["RU"];
            if (args.Count() == 2)
            {
                Console.WriteLine("*****reading the input from disk*****");
                Console.WriteLine("reading \"" + args[0] + "\"........");
                fileDetach=args[1];
            }
                //Console.WriteLine("reading \"" + args[1] + "\"........");
                inputDataShort = DataIO.ReadDataTable(fileDetach , true, '\t', 0);
                shortT_D = inputDataShort["time"];
                shortR_D = inputDataShort["RU"];
                Console.WriteLine("Done");
            //}
            SprModel sprm=new SprModel(shortT_A,shortR_A, shortT_D, shortR_D);
            sprm.SetAllParamters(_ka, _kd, _kM, _conc, _Rmax, _R0, var);

            List<int> lstFunc = new List<int>();
            lstFunc.Add(0/*ka*/); 
            lstFunc.Add(1/*kd*/);
            //lstFunc.Add(2/*kM*/);
            lstFunc.Add(3/*conc*/);
            lstFunc.Add(4/*Rmax*/);
            lstFunc.Add(5);/*R0*/
            lstFunc.Add(6/*var*/);
            sprm.setFunctionDelegateForUpdating(lstFunc);

            Console.WriteLine("setting up the lower/upper bounds of parameters..........");
            List<List<double>> bounds = new List<List<double>>();
            bounds.Add(new List<double> { 1E3, 1E10 });//ka, (1, 1E16)
            bounds.Add(new List<double> { 0, 10E0 });//kd, (0, 100)
            //bounds.Add(new List<double> { 1E5, 1E11 });//kM, (0, 1E15)
            bounds.Add(new List<double> { 0, 1E-4 });//conc, (0, 10)
            bounds.Add(new List<double> { 5E0, 10E1 });//Rmax, (0, 1E5)
            bounds.Add(new List<double> { 1, 10E1 });//R0, (0, 1E5)
            bounds.Add(new List<double> { 0, 3});//var, (0, 3000)

            Console.WriteLine("building the gibbs sampler...........");
            GibbsSampler.GibbsSampler gb = new GibbsSampler.GibbsSampler(/*init values*/new List<double> { 1E6, 1E-4, /*8E9,*/1E-8, 33.3462472360981,33.949870936, 1.5541821 /* */}, sprm.updateFunctionDistribution, bounds);
            //GibbsSampler.GibbsSampler gb = new GibbsSampler.GibbsSampler(new List<double> {1545903.57844479,	0.000322900628507969,	9.91536866021184E-09,	29.8920610679908,	29.5843847101981,	1.65389525923268}, sprm.updateFunctionDistribution, bounds);
            Console.WriteLine("now ready to run estimation........");
            List<List<double>> output = new List<List<double>>();
            output = gb.Run(15000);


            Console.WriteLine("writing the output file.......");
            StreamWriter writer = new StreamWriter("learReg.txt");
            //writer.WriteLine("ka\tkb\tkM\tconc\tRmax\tR0\tVar");
            for (int i = 0; i < output[0].Count; i++)
            {
                for (int j = 0; j < output.Count; j++)
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

            Console.WriteLine("Done...........");

        }//end of main

        static void writeOutput(List<double> lst1, List<double> lst2, string _filename, List<string> header)
        {
            StreamWriter writer = new StreamWriter(_filename);
            writer.WriteLine(header[0] + "\t" + header[1]);
            for (int i = 0; i < lst1.Count; i++)
            {
                writer.WriteLine(lst1[i] + "\t" + lst2[i]);
            }
            writer.Close();

        }

        //define variables.
        private static NumericalIntegrationOfDynamics nid;
        private static double _ka;
        private static double _kd;
        private static double _conc;
        private static double _Rmax;
        private static double _duration_attach;
        private static double _duration_detach;
        private static double _kM;
        private static double _deltaT;
        private static double _R0;
        
    }//end of class
}//end of namespace
