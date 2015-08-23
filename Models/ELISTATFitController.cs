using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BayesianEstimateLib;
namespace Models
{
    public class ELISTATFitController:FitController 
    {
        public ELISTATFitController()
        {
            //everything is null so far
        }

        public ELISTATFitController(List<List<double>> _X, List<double> _Y ):base(_X, _Y)
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
            C_Model = new ELISTATModel(new List<double> { 0.00001, 2.0, 0.0001, 0.1 });
            this.Read("ELISTAT_2014116_Abs_BAP0105.txt");
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
            C_Model.SetAllXs(this.C_X );
            C_Model.SetAllYs(this.C_Y);

            //set up prior
            List<List<double>> prior = new List<List<double>> (4);
            prior.Add(new List<double>() { 1e-5 ,1E4 });
            for(int i=0;i<2;i++)
            {
                prior.Add(new List<double>(){0,1E4});
                
            }
            prior.Add(new List<double>() { 0.01, 0.01 });

            C_Model.SetupPrior(prior);
            //set up bounds
            List<List<double>> bounds = new List<List<double>>(4);
            bounds.Add(new List<double>() { 0, 1E-7 });
            for (int i = 0; i < 2; i++)
            {
                bounds.Add(new List<double>() { 0, 4 });

            }
            bounds.Add(new List<double>() { 0, 1 });

            C_Model.SetupParameterBounds(bounds);
            this.C_Bounds = bounds;

            //set up parameter initials
            this.C_Parameters = new List<double> { 1E-10,2.0, 0.0001, 0.0001 };

            //set up parameter for updating list
            List<int> lstFunc = new List<int>();
            lstFunc.Add(0 /*KD*/);
            lstFunc.Add(1 /*a*/);
            lstFunc.Add(2/*b*/);
            lstFunc.Add(3/*var*/);
            C_Model.setFunctionDelegateForUpdating(lstFunc);
        }
        /// <summary>
        /// not implemented so far
        /// </summary>
        public override void Read(string _fileName)
        {
            Console.WriteLine("Start reading the file.........");
            Dictionary<int, List<double>> dt=DataIO.ReadDataTable(_fileName);
            List<double> temp = dt[0];
            C_X = new List<List<double>>();
            for (int i=0; i < temp.Count; i++)
            {
                C_X.Add(new List<double>() { temp[i] });
            }
            C_Y = dt[2];
        }
    }
}
