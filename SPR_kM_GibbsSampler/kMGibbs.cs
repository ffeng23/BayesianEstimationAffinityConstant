using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using NelderMeadMethod;

using AccessoryLib;

namespace SPR_kM_GibbsSampler
{
    /// <summary>
    /// caller to run Gibbs sampler to estimate the kM/Ckc from the training dataset
    /// </summary>
    class kMGibbs
    {
        static void Main(string[] args)
        {
            int totalNumOfInputs = 1; //there is one set of input files, need to be set up correctly if there is more than one set
            
            Console.WriteLine("run simulation....................");
            //setting up the input for the parameters
            _ka = new List<double>(totalNumOfInputs ) ;
            _ka.Add(1.1E5);//this has to be set up correctly too.

            _kd = new List<double>(totalNumOfInputs) ;
            _kd.Add( 1.28E-4);

            _conc = new List<double>(totalNumOfInputs);
            _conc.Add( 4E-8);

            _Rmax = new List<double>(totalNumOfInputs);
            _Rmax.Add( 152);

            _R0=new List<double>(totalNumOfInputs);
            _R0.Add( 97.15);
            //_duration_attach = 400;
            //_duration_detach = 400;
            C_var = new List<double>(totalNumOfInputs);
            C_var.Add( 2);
            
            
            C_Ckc = 1E15; //this is now the initial values
            C_flowRate = new List<double>(totalNumOfInputs);
            C_flowRate.Add( 50);

            C_molecularWeight = 15.6E3;
            
            Console.WriteLine("Reading sensorgram data....................");

            //now we need to shorten the data input
            List<List<double>> shortT_A = new List<List<double>>(totalNumOfInputs );
            List<List<double>> shortR_A = new List<List<double>>(totalNumOfInputs);
            List<List<double>> shortT_D = new List<List<double>>(totalNumOfInputs);
            List<List<double>> shortR_D = new List<List<double>>(totalNumOfInputs);
            
            //testing read in data
            Dictionary<string, List<double> > inputDataShort;//=BayesianEstimateLib.DataIO.ReadDataTable("simulation_attach_noiseShort.txt", true, '\t', 0);


            //***********************************************GIBBS SAMPLER********************
            //now testing Gibbs sampler


            //look to see whether we want to use the new simulated input or the one from the file
            List<string> fileAttach =new List<string>(totalNumOfInputs);
            fileAttach.Add( "attachingData0_240secL1A1AnalyzedIL-2Demo.txt");
            List<string> fileDetach =new List<string>(totalNumOfInputs);
            fileDetach.Add("detachingData0_240secL1A1AnalyzedIL-2Demo.txt");
            if (args.Count() == 2)
            {
                Console.WriteLine("*****reading the input file name from Command line*****");
                
                fileAttach[0]=args[0];
                fileDetach[0]=args[1];
            }
            for (int i = 0; i < totalNumOfInputs; i++)
            {
            
                //we need to read in the input
                Console.WriteLine("reading \"" + fileAttach[i] + "\"........");
                inputDataShort = DataIO.ReadDataTable(fileAttach[i], true, '\t', 0);
                shortT_A.Add(inputDataShort["time"]);
                shortR_A.Add(inputDataShort["RU"]);
            
                Console.WriteLine("reading \"" + fileDetach[i] + "\"........");
                inputDataShort = DataIO.ReadDataTable(fileDetach[i] , true, '\t', 0);
                shortT_D.Add(inputDataShort["time"]);
                shortR_D.Add( inputDataShort["RU"]);
                
            }
            Console.WriteLine("Done");
            
            //building the models to run the Gibbs sampler
            Console.WriteLine("Building models to run Gibbs sampler.......");
            kMEstimationModel  sprm=new kMEstimationModel (shortT_A,shortR_A, shortT_D, shortR_D);
            sprm.SetkMRelatedParameters(C_molecularWeight, C_flowRate, 25);
            
            sprm.SetAllParamters(_ka, _kd,  _conc, _Rmax, _R0, C_var,C_Ckc);

            List<int> lstFunc = new List<int>();
            lstFunc.Add(0/*Ckc*/);
            lstFunc.Add(1/*Rmax*/);
            lstFunc.Add(2/*R0*/);
            lstFunc.Add(3/*var*/);
            sprm.setFunctionDelegateForUpdating(lstFunc);

            Console.WriteLine("setting up the lower/upper bounds of parameters..........");
            List<List<double>> bounds = new List<List<double>>(); //one set of 
            bounds.Add(new List<double> { 1E10, 1E21 });//Ckc, (1, 1E16)
            bounds.Add(new List<double> { 0, 5E2 });//Rmax, (0, 100)
            bounds.Add(new List<double> { 0, 5E2 });//R0, (0, 1E15)
            bounds.Add(new List<double> { 0, 1E1 });//conc, (0, 10)
            //bounds.Add(new List<double> { 5E0, 5E2 });//Rmax, (0, 1E5)
            //bounds.Add(new List<double> { 1, 1E3 });//R0, (0, 1E5)
            //bounds.Add(new List<double> { 0, 20 });//var, (0, 3000)

            Console.WriteLine("building the gibbs sampler...........");
            List<double> initials = new List<double>();
            initials.Add(1E15);
            initials.Add(160);
            initials.Add(100);
            initials.Add(3);
            GibbsSampler.GibbsSampler gb = new GibbsSampler.GibbsSampler(initials, sprm.updateFunctionDistribution, bounds);
            //GibbsSampler.GibbsSampler gb = new GibbsSampler.GibbsSampler(new List<double> {1545903.57844479,	0.000322900628507969,	9.91536866021184E-09,	29.8920610679908,	29.5843847101981,	1.65389525923268}, sprm.updateFunctionDistribution, bounds);
            Console.WriteLine("now ready to run estimation........");
            List<List<double>> output = new List<List<double>>();
            output = gb.Run(5000);


            Console.WriteLine("writing the output file.......");
            StreamWriter writer = new StreamWriter("learReg.txt");
            //writer.WriteLine("Ckc\tVar");
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
        //private static NumericalIntegrationOfDynamics nid;
        private static List<double> _ka;
        private static List<double> _kd;
        private static List<double> _conc;
        private static List<double> _Rmax;
        private static List<double> _duration_attach;
        private static List<double> _duration_detach;
        private static List<double> _kM;
        //private static double _deltaT;
        private static List<double> _R0;
        private static List<double> C_var;

        private static double C_Ckc;
        private static double C_molecularWeight;
        private static List<double> C_flowRate;
        
    }//class
}//namespace
