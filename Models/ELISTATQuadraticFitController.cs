using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BayesianEstimateLib;

using AccessoryLib;

namespace Models
{
    public class ELISTATQuadraticFitController:FitController
    {
        public ELISTATQuadraticFitController()
        {
            //everything is null so far
        }

        public ELISTATQuadraticFitController(List<List<double>> _X, List<double> _Y ):base(_X, _Y)
        {
            //the following code was done in the base class
            /*
            this.C_Model = null;
            this.C_X=_X;
            this.C_Y=_Y;

            C_Parameters=null;
            C_Bounds = null;
             */
        }
        /// <summary>
        /// set up parameters, bounds and prior information 
        /// </summary>
        public override void SetupModel()
        {
            //need to set up parameter
            C_Model = new ELISTATQuadraticModel(new List<double> { 1.9e-11, 3e-21, 1e21, 0.00, 0.0001 });
            this.Read("simulatedDataAg3E-9Ab60E-9KD1.3E-9-totalNoise.txt");
            //this.Read("ELISTAT_format_plateA_nov13Run.txt");
            //this.Read("simulatedDataAg3E-9Ab5E-9KD1.3E-9-totalNoise.txt");
            //this.Read("simulatedDataAg3E-9Ab0.3E-9KD1.3E-9-totalNoise.txt");
            //this.Read("simulatedDataAg36E-9Ab0.36E-9KD1.3E-9-total.txt");
            //this.Read("ELISTAT_format_plateA_nov13Run1.txt");
            /*List<List<double>> Xsim = new List<List<double>>(200);
            for (int i = 0; i < 200; i++)
            {
                for (int j = 0; j < 1; j++)
                {
                    Xsim.Add( new List<double> { -100+ i });
                }
            }

            List<double> Ysim=C_Model.SimulateYValues(Xsim, 1.5);
            */
            C_Model.SetAllXs(this.C_X);
            C_Model.SetAllYs(this.C_Y);

            //set up prior
            List<List<double>> prior = new List<List<double>>(5);
            prior.Add(new List<double>() { 1e-11, 1E5 });
            prior.Add(new List<double>() { 1e-15, 1E10 });

            prior.Add(new List<double>() { 1e15, 1E10 });
            prior.Add(new List<double>() { 0, 1E4 });
            prior.Add(new List<double>() { 0.01, 0.01 });

            C_Model.SetupPrior(prior);
            //set up bounds
            List<List<double>> bounds = new List<List<double>>();
            bounds.Add(new List<double>() { 1E-13, 1E-6 });
            bounds.Add(new List<double>() { 1E-20, 1E-6 });
            
                bounds.Add(new List<double>() { 1E4, 1E20 });

                //bounds.Add(new List<double>() { 0, 1E3 });

            bounds.Add(new List<double>() { 0, 1 });

            C_Model.SetupParameterBounds(bounds);
            this.C_Bounds = bounds;

            //set up parameter initials
            this.C_Parameters = new List<double> { 1.3e-11, 3e-15, 1e15,/* 0.00000,*/ 0.001 };

            //set up parameter for updating list
            List<int> lstFunc = new List<int>();
            lstFunc.Add(0 /*KD*/);
            lstFunc.Add(1 /*Ag_t*/);
            lstFunc.Add(2/*a*/);
            //lstFunc.Add(3/*b*/);
            lstFunc.Add(4/*var*/);
            C_Model.setFunctionDelegateForUpdating(lstFunc);
        }
        /// <summary>
        /// not implemented so far
        /// </summary>
        public override void Read(string _fileName)
        {
            Console.WriteLine("Start reading the file.........");
            Dictionary<int, List<double>> dt = DataIO.ReadDataTable(_fileName);
            List<double> temp = dt[0];
            C_X = new List<List<double>>();
            for (int i = 0; i < temp.Count; i++)
            {
                C_X.Add(new List<double>() { temp[i] });
            }
            C_Y = dt[2];
        }

    }//end of class
}
