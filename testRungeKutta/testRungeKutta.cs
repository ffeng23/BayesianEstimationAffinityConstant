﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using RungeKuttaMethod;
using System.IO;
using BayesianEstimateLib;
namespace testRungeKutta
{
    class testRungeKutta
    {
        static void Main(string[] args)
        {
            //RungeKutta rk = new RungeKutta();
            int count=1000;
            List<double> input=new List<double>(count);
            double step=1;
            for(int i=0;i<count;i++)
            {
                input.Add(i*step);
            }
            List<double> output=RungeKutta.Solution(Quadratic_Derivatives, new List<double>{input[0], input[input.Count-1]}, step, 2.0);

            BayesianEstimateLib.DataIO.WriteDataTable(input, output,"rungeKutta_der.txt",new List<string>{"x", "y"});

            List<double> outputOrig = new List<double>();
            for (int i = 0; i < count; i++)
            {
                outputOrig.Add(Quadratic(input[i], input[i]));
            }
            BayesianEstimateLib.DataIO.WriteDataTable(input, outputOrig, "rungeKutta_ori.txt", new List<string> { "x", "y" });

            List<double> eul = Euler(new List<double> { input[0], input[input.Count - 1] }, 2.0, step);
            BayesianEstimateLib.DataIO.WriteDataTable(input, eul, "Euler.txt", new List<string> { "x", "y" });

            //***********testing new function and with the two different implementation of RK4
            count = 100;
            List<double> inputNew = new List<double>(count);
            step = 0.05;
            for (int i = 0; i < count; i++)
            {
                inputNew.Add(i * step);
            }
            List<double> outputNew = RungeKutta.Solution(Function_Derivatives, new List<double> { input[0], input[input.Count - 1] }, step, 0.5);

            BayesianEstimateLib.DataIO.WriteDataTable(inputNew, outputNew, "rungeKutta_derNew.txt", new List<string> { "x", "y" });

            outputNew = RungeKutta.SolutionH(Function_Derivatives, new List<double> { input[0], input[input.Count - 1] }, step, 0.5);

            BayesianEstimateLib.DataIO.WriteDataTable(inputNew, outputNew, "rungeKutta_derNewH.txt", new List<string> { "x", "y" });

            outputOrig = new List<double>();
            for (int i = 0; i < count; i++)
            {
                outputOrig.Add(functionOrig (inputNew[i]));
            }
            BayesianEstimateLib.DataIO.WriteDataTable(inputNew, outputOrig, "rungeKutta_oriNew.txt", new List<string> { "x", "y" });

            //*************
            //now we are testing the spr data, and compare the result between the runge kutta and the langmuir model ones
            double _ka = 1E5;
            double _kd = 1E-3;
            double _conc = 3E-8;
            double _Rmax = 31.5;

            double _duration_attach = 400;
            double _duration_detach = 400;

            double _kM = 3.15E8; //this is now in RU/S
            double _deltaT = 0.001;
            double _R0 = 20;

            Console.WriteLine("run simulation....................");

            NumericalIntegrationOfDynamics nid = new NumericalIntegrationOfDynamics(_ka, _kd, _conc, _Rmax, _R0, _duration_attach, _duration_detach, _kM, _deltaT);
            
            List<string> header = new List<string>(2);
            header.Add("time");
            header.Add("RU");
            //doing the simulation for different Rmax

            nid.run_Attach();
            nid.R0 = _R0;
            nid.run_Detach();


            Console.WriteLine("write output..................");
            DataIO.WriteDataTable(nid.Time_Attach, nid.RU_Attach, "simulation_attachRK.txt", header);
            DataIO.WriteDataTable(nid.Time_Detach, nid.RU_Detach, "simulation_detachRK.txt", header);

            nid.run_AttachEuler();
            nid.R0 = _R0;
            nid.run_DetachEuler();
            Console.WriteLine("write output..................");
            DataIO.WriteDataTable(nid.Time_Attach, nid.RU_Attach, "simulation_attachEU.txt", header);
            DataIO.WriteDataTable(nid.Time_Detach, nid.RU_Detach, "simulation_detachEU.txt", header);
            Console.WriteLine("Done..................");
            Langmuir lid = new Langmuir(_ka, _kd, _conc, _Rmax, _R0, _duration_attach, _duration_detach, _deltaT);
            lid.run_Attach();
            lid.R0 = _R0;
            lid.run_Detach();
            Console.WriteLine("write output..................");
            DataIO.WriteDataTable(lid.Time_Attach, lid.RU_Attach, "simulation_attachLG.txt", header);
            DataIO.WriteDataTable(lid.Time_Detach, lid.RU_Detach, "simulation_detachLG.txt", header);
            Console.WriteLine("Done..................");
        }
        static double Quadratic_Derivatives(double _x, double _y)
        {
            //y=(x-a)^2+b
            return 2 * _x;
        }
        static double Function_Derivatives(double _t, double _y)
        {
            //y=(x-a)^2+b
            return _y-_t* _t+1;
        }
        static double Quadratic(double _x, double _y)
        {
            return _x * _x + 2.0;
        }
        static double functionOrig(double _x)
        {
            return _x * _x + 2.0*_x+1-0.5*Math.Exp(_x);
        }

        static List<double> Euler(List<double> _input, double _init, double _h)
        {
            List<double> ret = new List<double>();
            ret.Add(_init);
            for (int i = 1; i <= (_input[1] - _input[0]) / _h; i++)
            {
                double temp=ret[i-1]+_h*Quadratic_Derivatives(_h*i+_input[0],_h*i+_input[0]);
                ret.Add(temp);
            }

            return ret;
        }
    }

    
}
