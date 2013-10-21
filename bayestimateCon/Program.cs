using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BayesianEstimateLib;
using Savitzky_GolaySmoothingLib;

using System.IO;

namespace bayestimateCon
{
    class BayEstimateCon
    {
        static void Main(string[] args)
        {
            int MC_TotalCycles = 100000;
            //Console.WriteLine("args.count " + args.Count() + " args[0] " + args[0]);
            if(args.Count()>=1)
            {
                MC_TotalCycles = Convert.ToInt32(args[0]);
            }
            //getting values first
            _ka = 6E6;
            _kd = 1E-3;
            _conc = 3E-8;
            _Rmax = 31.5;
            
            _duration_attach = 400;
            _duration_detach = 400;
            
            _kM = 3.15E7; //this is now in RU/S
            _deltaT = 0.1;

            Console.WriteLine("run simulation....................");

            nid = new NumericalIntegrationOfDynamics(_ka, _kd, _conc, _Rmax,200,_duration_attach, _duration_detach,_kM, _deltaT);
            double r0=0;
            List<string> header = new List<string>(2);
            header.Add("time");
            header.Add("RU");
            //doing the simulation for different Rmax
            int counts = 1; 
            
            //int frameSize = 30, degreeOfPolynomial = 2;
            //SG_Smoothing sgs = new SG_Smoothing(frameSize, degreeOfPolynomial);
            for (int i = 0; i < counts; i++)
            {
                nid.run_Attach();
                nid.R0 = nid.RU_Attach[nid.RU_Attach.Count - 1];
                nid.run_Detach();
                

                Console.WriteLine("write output..................");
                writeOutput(nid.Time_Attach, nid.RU_Attach, "simulation_attach.txt", header);
                writeOutput(nid.Time_Detach, nid.RU_Detach, "simulation_detach.txt", header);


                //writeOutput(nid.Time_Attach, nid.RU_Attach, "simulation_attach_Rmax"+(_Rmax-stepOfRmax*i) +".txt", header);
                //writeOutput(nid.Time_Detach, nid.RU_Detach, "simulation_detach_Rmax" + (_Rmax - stepOfRmax*i) + ".txt", header);

                Console.WriteLine("add noise...............");
                if(i==0)
                    r0 = nid.RU_Attach[nid.RU_Attach.Count - 1];
                nid.addNoise(1.5);
                writeOutput(nid.Time_Attach, nid.RU_Attach, "simulation_attach_noise.txt", header);
                writeOutput(nid.Time_Detach, nid.RU_Detach, "simulation_detach_noise.txt", header);
                
                //writeOutput(nid.Time_Attach, nid.RU_Attach, "simulation_attach_noise_Rmax" + (_Rmax - stepOfRmax*i) + ".txt", header);
                //writeOutput(nid.Time_Detach, nid.RU_Detach, "simulation_detach_noise_Rmax" + (_Rmax - stepOfRmax*i) + ".txt", header);
                
                //get short list of attach and detach
                /*Dictionary<string, List<Double>> dt_short = new Dictionary<string, List<double>>();
                //List<int> key = dt.Keys.ToList();
                dt_short.Add("Time", new List<double>());
                dt_short.Add("RU", new List<double>());
                for (int j = 0; j < nid.Time_Attach.Count; j = j + 900)
                {
                    dt_short["Time"].Add(nid.Time_Attach[j]);
                    dt_short["RU"].Add(nid.RU_Attach[j]);
                }


                //now doing the smoothing, hopefully this will reduce the 
                List<double> dt_smoothed = sgs.Smooth(0, dt_short["RU"]);
                List<string> header1 = new List<string> { "time", "RU" };
                DataIO.WriteDataTable(dt_short["Time"] , dt_smoothed, "simulation_attach_noise_smoothed_Rmax" + (_Rmax - stepOfRmax*i) + ".txt", header1);

                //now doing slope estimation
                List<double> dt_slope = sgs.Smooth(1, dt_short["RU"]);
                header1 = new List<string> { "time", "RU" };
                DataIO.WriteDataTable(dt_short["Time"], dt_slope, "simulation_attach_noise_slope_Rmax" + (_Rmax - stepOfRmax * i) + ".txt", header1);

                //set the nid for next round
                nid.setParameters(_ka, _kd, _kM, _conc, _Rmax - stepOfRmax*(i+1), nid.R0);*/
            }
            //nid.setParameters(_ka, _kd, _kM, _conc, _Rmax, nid.R0);
            Console.WriteLine("done............");

            //using the noise array to do the estimation
            Console.WriteLine("doing MCMC estimation............");
            //LogTextBlock.Text += "\rdoing MCMC...........";
            MCMC mc = new MCMC_MH ();
            mc.getInput(nid.Time_Attach, nid.RU_Attach,nid.Time_Detach, nid.RU_Detach);
            mc.MC_TotalCycles = MC_TotalCycles;
            //mc.R00 = r0;
            //Console.WriteLine("R0 is " + r0);

            mc.run();
            Console.WriteLine("writing output.........");
            mc.writeOutput();
            //LogTextBlock.Text += "\rdone.....";
            Console.WriteLine("done..........");

            //Console.WriteLine("R0 is " + r0);
            //doing the simultion using the means
            _ka = mc.mean_ka ;
            _kd = mc.mean_kd;
            _conc = mc.mean_conc ;
            _Rmax = mc.mean_Rmax ;
            _R0 = mc.mean_r0;
            _duration_attach = 400;
            _duration_detach = 1300;
            _kM = mc.mean_kM ;
            nid = new NumericalIntegrationOfDynamics(_ka, _kd, _conc, _Rmax, _R0, _duration_attach, _duration_detach,  _kM, _deltaT);
            nid.run_Attach();

            Console.WriteLine("write output for final simultion..................");

            writeOutput(nid.Time_Attach, nid.RU_Attach, "simulationMean.txt", header);
            nid.run_Detach();
            writeOutput(nid.Time_Detach, nid.RU_Detach, "simulationMean_detach.txt", header);

            Console.WriteLine("add noise...............");

//nid.addNoise(mc.mean_sigma);
            //writeOutput(nid.Time_Attach, nid.RU_Attach, "simulation_noiseMean.txt", header);
            Console.WriteLine("done............");


            //*************************
            //doing the simultion using the last one
            Console.WriteLine("now we are ready to do the simulation with last sets of parameters............");

            _ka = mc._ka;
            _kd = mc._kd;
            _conc = mc._conc;
            _Rmax = mc._Rmax;
            _R0 = mc.MCMC_r0;

            _duration_attach = 400;
            _duration_detach = 1300;
            
            _kM = mc._kM;
            Console.WriteLine("_conc<-" + _conc);
            Console.WriteLine("_ka<-" + _ka);
            Console.WriteLine("_kd<-" + _kd);
            Console.WriteLine("_kM<-" + _kM);
            Console.WriteLine("_sigma<-" + mc._sigma);
            Console.WriteLine("Rmax<-" + _Rmax);
            Console.WriteLine("R0<-" + _R0);


            nid = new NumericalIntegrationOfDynamics(_ka, _kd, _conc, _Rmax, _R0, _duration_attach, _duration_detach, _kM, _deltaT);
            nid.run_Attach();

            Console.WriteLine("write output for Last simultion..................");

            writeOutput(nid.Time_Attach, nid.RU_Attach, "simulationLast.txt", header);
            nid.run_Detach();
            writeOutput(nid.Time_Detach, nid.RU_Detach, "simulationLast_detach.txt", header);

            Console.WriteLine("write output for MCMC Last cycle simultion..................");

            writeOutput(mc.MC_nid.Time_Attach, mc.MC_nid.RU_Attach, "simulationMCLast.txt", header);
            
            //for simulatoin of last next parameters
            nid = new NumericalIntegrationOfDynamics(mc.nextParameters["ka"], mc.nextParameters["kd"], mc.nextParameters["conc"], mc.nextParameters["Rmax"],mc.nextParameters["R0"], 
                _duration_attach, _duration_detach,  mc.nextParameters["kM"], _deltaT);
            nid.run_Attach();

            Console.WriteLine("write output for Last next simultion..................");

            writeOutput(nid.Time_Attach, nid.RU_Attach, "simulationLastNext.txt", header);


            Console.WriteLine("_conc<-" + mc.nextParameters["conc"]);
            Console.WriteLine("_ka<-" + mc.nextParameters["ka"]);
            Console.WriteLine("_kd<-" + mc.nextParameters["kd"]);
            Console.WriteLine("_kM<-" + mc.nextParameters["kM"]);
            Console.WriteLine("_sigma<-" + mc.nextParameters["sigma"]);
            Console.WriteLine("Rmax<-" + mc.nextParameters["Rmax"]);
            Console.WriteLine("R0<-" + mc.nextParameters["R0"]); 
            //nid.addNoise(mc._sigma);
            //writeOutput(nid.Time_Attach, nid.RU_Attach, "simulation_noiseMean.txt", header);
            Console.WriteLine("done............");


        }

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
        
        private List<double> _time;
        private List<double> _ru;

    }//end of class.
}
