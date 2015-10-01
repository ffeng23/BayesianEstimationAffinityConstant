using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BayesianEstimateLib;
using AccessoryLib;
namespace Models
{
    /// <summary>
    /// this is the fit controller to call to do the simple quadratic equation fitting for testing purposes.
    /// </summary>
    public class QuadraticFitController:FitController
    {
        public QuadraticFitController()
        {
            //everything is null so far
        }

        public QuadraticFitController(List<List<double>> _X, List<double> _Y ):base(_X, _Y)
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
            C_Model = new QuadraticModel(new List<double> { 3000.0, 1500, 1.5 });
            List<List<double>> Xsim = new List<List<double>>(200);
            for (int i = 0; i < 200; i++)
            {
                for (int j = 0; j < 1; j++)
                {
                    Xsim.Add( new List<double> { -100+ i });
                }
            }

            List<double> Ysim=C_Model.SimulateYValues(Xsim, 1.5);
            C_Model.SetAllXs(Xsim);
            C_Model.SetAllYs(Ysim);

            //write to the disk the Observed data
            DataIO.WriteDataTable(Ysim, Xsim, "QuadraticData.txt", new List<string> { "Y", "X" });
            //set up prior
            List<List<double>> prior = new List<List<double>> (3);
            for(int i=0;i<2;i++)
            {
                prior.Add(new List<double>(){0,1E6});
                
            }
            prior.Add(new List<double>() { 0.01, 0.01 });

            C_Model.SetupPrior(prior);
            //set up bounds
            List<List<double>> bounds = new List<List<double>>(3);
            for (int i = 0; i < 2; i++)
            {
                bounds.Add(new List<double>() { -1E4, 1E4 });

            }
            bounds.Add(new List<double>() { 0, 100 });

            C_Model.SetupParameterBounds(bounds);
            this.C_Bounds = bounds;

            //set up parameter initials
            this.C_Parameters = new List<double> { 300, 150, 10 };

            //set up parameter for updating list
            List<int> lstFunc = new List<int>();
            lstFunc.Add(0 /*a*/);
            lstFunc.Add(1 /*b*/);
            lstFunc.Add(2/*var*/);
            C_Model.setFunctionDelegateForUpdating(lstFunc);
        }
        /// <summary>
        /// not implemented so far
        /// </summary>
        public override void Read(string _fileName)
        {
            Console.WriteLine("We don't implement in this module, return........");
        }
        /*run was called in the base class
        /// <summary>
        /// the work horse to run Gibbs sampler 
        /// </summary>
        /// <param name="_NumSteps">number of steps for sampler to draw</param>
        /// <returns></returns>
        //public List<List<double>> Run(int _NumSteps)
        {
            if(C_Parameters!=null)
            {
                Console.WriteLine("********ERROR*******:the parameters have not been set up, try again........");
                return null;
            }
            GibbsSampler.GibbsSampler gbs = new GibbsSampler.GibbsSampler(C_Parameters, C_Model.updateFunctionDistribution, C_Bounds);
            return gbs.Run(_NumSteps );
        }
        */

    }//end of class
}
