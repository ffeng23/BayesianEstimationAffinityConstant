using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using GibbsSampler;

namespace Models
{
    /// <summary>
    /// this is template class to direct the fitting, basically we need to set up the model
    /// </summary>
    public abstract class FitController
    {
        public FitController()
        {
            this.C_Model = null;
            this.C_X = null;
            this.C_Y = null;

            C_Parameters = null;
            C_Bounds = null;
        }
        public FitController(List<List<double>> _X, List<double> _Y )
        {
            this.C_Model = null;
            this.C_X=_X;
            this.C_Y=_Y;

            C_Parameters=null;
            C_Bounds = null;
        }
        /// <summary>
        /// set up parameters, bounds and prior information 
        /// </summary>
        public abstract void SetupModel();

        /// <summary>
        /// the work horse to run Gibbs sampler 
        /// </summary>
        /// <param name="_NumSteps">number of steps for sampler to draw</param>
        /// <returns></returns>
        public List<List<double>> Run(int _NumSteps)
        {
            if(C_Parameters==null||C_Bounds==null||C_Model.ParameterPrior==null)
            {
                Console.WriteLine("********ERROR*******:the parameters have not been set up, try again........");
                return null;
            }

            GibbsSampler.GibbsSampler gbs = new GibbsSampler.GibbsSampler(C_Parameters, C_Model.updateFunctionDistribution, C_Bounds);
            return gbs.Run(_NumSteps );
        }

        public abstract void Read(string _fileName);

        protected Model C_Model; //the mathematic model to fit
        protected List<double> C_Parameters;
        protected List<List<double>> C_Bounds; //2-D array holding the bounds.
        //protected GibbsSampler.GibbsSampler C_GS;
        protected List<List<double>> C_X;
        protected List<double> C_Y;


    }
}
