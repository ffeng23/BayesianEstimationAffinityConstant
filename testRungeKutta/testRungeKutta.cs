using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using RungeKuttaMethod;
using System.IO;
using BayesianEstimateLib;
using AccessoryLib;
namespace testRungeKutta
{
    class testRungeKutta
    {
        static void Main(string[] args)
        {

            Console.WriteLine("Start doing test........");
            //RungeKutta rk = new RungeKutta();
            int count=1000;
            List<double> input=new List<double>(count);
            double step=1;
            for(int i=0;i<count;i++)
            {
                input.Add(i*step);
            }
            List<double> output=RungeKutta.Solution(Quadratic_Derivatives, new List<double>{input[0], input[input.Count-1]}, step, 2.0);

            DataIO.WriteDataTable(input, output,"rungeKutta_der.txt",new List<string>{"x", "y"});

            List<double> outputOrig = new List<double>();
            for (int i = 0; i < count; i++)
            {
                outputOrig.Add(Quadratic(input[i], input[i]));
            }
            DataIO.WriteDataTable(input, outputOrig, "rungeKutta_ori.txt", new List<string> { "x", "y" });

            List<double> eul = Euler(new List<double> { input[0], input[input.Count - 1] }, 2.0, step);
            DataIO.WriteDataTable(input, eul, "Euler.txt", new List<string> { "x", "y" });

            //***********testing new function and with the two different implementation of RK4
            count = 100;
            List<double> inputNew = new List<double>(count);
            step = 0.05;
            for (int i = 0; i < count; i++)
            {
                inputNew.Add(i * step);
            }
            List<double> outputNew = RungeKutta.Solution(Function_Derivatives, new List<double> { input[0], input[input.Count - 1] }, step, 0.5);

            DataIO.WriteDataTable(inputNew, outputNew, "rungeKutta_derNew.txt", new List<string> { "x", "y" });

            outputNew = RungeKutta.SolutionH(Function_Derivatives, new List<double> { input[0], input[input.Count - 1] }, step, 0.5);

            DataIO.WriteDataTable(inputNew, outputNew, "rungeKutta_derNewH.txt", new List<string> { "x", "y" });

            outputOrig = new List<double>();
            for (int i = 0; i < count; i++)
            {
                outputOrig.Add(functionOrig (inputNew[i]));
            }
            DataIO.WriteDataTable(inputNew, outputOrig, "rungeKutta_oriNew.txt", new List<string> { "x", "y" });

            //*************
            //now we are testing the spr data, and compare the result between the runge kutta and the langmuir model ones
            double _ka = 1E5;
            double _kd = 1E-3;
            double _conc = 3E-8;
            double _Rmax = 31.5;

            double _duration_attach = 400;
            double _duration_detach = 400;

            double _kM = 3.15E8; //this is now in RU/S
            double _deltaT = 0.1;
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

            //testing two state SPR model with Euler and RK4
            InducedFit ts = new InducedFit();
            _ka = 1E6; _kd = 0.005;
            double _ka2 = 0.003, _kd2=0.002;
            _conc = 1E-7; _Rmax = 42;
            ts.setParameters(new double[] { _ka, _kd, _ka2, _kd2, _conc, _Rmax,-1,-1,1000,1000,1 });
            Console.WriteLine("Doing two state with Euler scheme..............");
            ts.run_Attach_RK();
            ts.run_Detach_RK();
            DataIO.WriteDataTable(ts.Time_Attach, ts.RU_Attach, "simulation_TwoState_attach_Euler.txt", header);
            DataIO.WriteDataTable(ts.Time_Attach, ts.RU_Detach, "simulation_TwoState_detach_Euler.txt", header);

            //ts.setParameters(new double[] { _ka, _kd, _ka2, _kd2, _conc, _Rmax, -1, -1, 1000, 1000, 0.01 });
            //testing conformational selection SPR model with Euler and RK4
            ConformationalSelection  cs = new ConformationalSelection();
            _ka = 0.002; _kd = 0.001;
            _ka2 = 1E4; _kd2 = 0.001;
            _conc = 1E-7; _Rmax = 42;
            cs.setParameters(new double[] { _ka, _kd, _ka2, _kd2, _conc, _Rmax, -1,-1, 1000, 1000, 0.1 });
            Console.WriteLine("Doing cs with Euler scheme..............");
            cs.run_Attach_RK();
            cs.run_Detach_RK();
            DataIO.WriteDataTable(cs.Time_Attach, cs.RU_Attach, "simulation_CS_attach_Euler.txt", header);
            DataIO.WriteDataTable(cs.Time_Attach, cs.RU_Detach, "simulation_CS_detach_Euler.txt", header);

            //ts.setParameters(new double[] { _ka, _kd, _ka2, _kd2, _conc, _Rmax, -1, -1, 1000, 1000, 0.01 });
            //testing conformational selection SPR model with Euler and RK4
            TwoStates tws = new TwoStates();
            _ka = 1e5; _kd = 0.005;
            double _ka_if = 0.03; double _kd_if = 0.002;
            double _kon_cs = 1E4; double _koff_cs= 0.001;
            double _ka_cs = 0.002; double _kd_cs = 0.001;
            _conc = 1E-7; _Rmax = 250;
            tws.setParameters(new double[] { _ka, _kd, _ka_if, _kd_if,_kon_cs, _koff_cs, _ka_cs, _kd_cs, _conc, _Rmax, -1, -1, -1,-1, 1000, 1000, 0.1 });
            Console.WriteLine("Doing two states with Euler scheme..............");
            tws.run_Attach_RK();
            tws.run_Detach_RK();
            DataIO.WriteDataTable(tws.Time_Attach, tws.RU_Attach, "simulation_TWS_attach_Euler.txt", header);
            DataIO.WriteDataTable(tws.Time_Attach, tws.RU_Detach, "simulation_TWS_detach_Euler.txt", header);

            Console.WriteLine("Done.............");
            Console.WriteLine("please type Enter to exit.........");
            Console.ReadLine();
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
