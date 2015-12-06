using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using BayesianEstimateLib;

using AccessoryLib;

namespace Models
{
    public class FitController_UnifiedTwoState:FitController 
    {
        public FitController_UnifiedTwoState(double _conc, double _attach_duration):base()
        {
            //everything is null so far

            //need to do more than base does
            this.C_Y_Detach = null;
            this.C_Conc = _conc;
            this.C_Duration_Attach = _attach_duration;
        }
        public void SetAttachDuration(double _attach_duration)
        {
            this.C_Duration_Attach = _attach_duration;
        }
        public void SetConc(double _conc)
        {
            this.C_Conc = _conc;
        }

        /*public ELISTATFitController(List<List<double>> _X, List<double> _Y ):base(_X, _Y)
        {
            //the following code was done in the base class
            / *
            this.C_Model = null;
            this.C_X=_X;
            this.C_Y=_Y;

            C_Parameters=null;
            C_Bounds = null;
             * /
        }*/
        /// <summary>
        /// set up parameters, bounds and prior information 
        /// </summary>
        public override void SetupModel()
        {
            this.Read("E:\\MSI_software\\ELIStat_HIP\\ELIStat\\SPR_twoState\\TwoStateInput_Full_noiseSD2.0.txt");
            List<List<double>> ys=new List<List<double>>();
            ys.Add(this.C_Y);
            ys.Add(this.C_Y_Detach);
            List<double> pms=new List<double> (){1e5 /*kon_if*/, 0.01 /*koff_if*/, 0.03 /*ka_if*/, 0.001 /*kd_if*/,
                    10000 /*kon_cs*/, 0.0008 /*koff_cs*/, 0.001 /*ka_cs*/, 0.001 /*kd_cs*/,
                    240 /*RMax*/, 2/*var*/
                };
            //need to set up parameter
            C_Model = new Model_UnifiedTwoState(pms,this.C_X, ys , 1E-7);
            
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
            //C_Model.SetAllXs(this.C_X );//set up the input Xs based on what is read in.
            //C_Model.SetAllYs(this.C_Y);//set up the input Ys based on what is read in.

            //set up prior
            List<List<double>> prior = new List<List<double>> (10);
            prior.Add(new List<double>() { 1e5, 1E8 });/*kon_if*/
            prior.Add(new List<double>() { 1e-4, 1E4 });/*koff_if*/
            prior.Add(new List<double>() { 1e-5, 1E4 });/*ka_if*/
            prior.Add(new List<double>() { 1e-5, 1E4 });/*kd_if*/

            prior.Add(new List<double>() { 1e4, 1E8 });/*kon_cs*/
            prior.Add(new List<double>() { 1e-5, 1E4 });/*koff_cs*/
            prior.Add(new List<double>() { 1e-5, 1E4 });/*ka_cs*/
            prior.Add(new List<double>() { 1e-5, 1E4 });/*kd_cs*/

            prior.Add(new List<double>() { 2e2, 1E4 });/*Rmax*/
            prior.Add(new List<double>() { 0.001, 0.001 });/*var*/
            /*
            for(int i=0;i<2;i++)
            {
                prior.Add(new List<double>(){0,1E4});
                
            }
            prior.Add(new List<double>() { 0.01, 0.01 });
            */
            C_Model.SetupPrior(prior);

            //set up bounds
            List<List<double>> bounds = new List<List<double>>(10);
            bounds.Add(new List<double>() { 0, 1E6 });/*kon_if*/
            bounds.Add(new List<double>() { 0, 1E-1 });/*koff_if*/
            bounds.Add(new List<double>() { 0, 1E-1 });/*ka_if*/
            bounds.Add(new List<double>() { 0, 1E-1 });/*kd_if*/

            bounds.Add(new List<double>() { 0, 1E6 });/*kon_cs*/
            bounds.Add(new List<double>() { 0, 1E-1 });/*koff_cs*/
            bounds.Add(new List<double>() { 0, 1E-1 });/*ka_cs*/
            bounds.Add(new List<double>() { 0, 1E-1 });/*kd_cs*/

            bounds.Add(new List<double>() { 0, 1E3 });/*Rmax*/
            bounds.Add(new List<double>() { 0, 3E3 });/*var*/
            
            C_Model.SetupParameterBounds(bounds);
            this.C_Bounds = bounds;

            //set up parameter initials
            this.C_Parameters = pms;// new List<double> { 1E-10, 2.0, 0.0001, 0.0001 };

            //set up parameter for updating list
            List<int> lstFunc = new List<int>();
            //lstFunc.Add(0 /*Kon_if*/);
            lstFunc.Add(1 /*Koff_if*/);
            lstFunc.Add(2/*Ka_if*/);
            lstFunc.Add(3/*kd_if*/);

            //lstFunc.Add(4 /*Kon_cs*/);
            lstFunc.Add(5 /*Koff_cs*/);
            lstFunc.Add(6/*Ka_cs*/);
            lstFunc.Add(7/*kd_cs*/);

            lstFunc.Add(8/*Rmax*/);
            lstFunc.Add(9/*var*/);

            C_Model.setFunctionDelegateForUpdating(lstFunc);

        }
        /// <summary>
        /// not implemented so far
        /// </summary>
        public override void Read(string _fileName)
        {
            Console.WriteLine("Start reading the file.........");
            Dictionary<int, List<double>> dt=DataIO.ReadDataTable(_fileName);
            //List<double> temp = dt[0];
            //List<double> tempY = dt[1];
            //now start parsing the attaching and detaching phases. 
            C_X = new List<List<double>>();
            List<double> tempX = dt[0];//first is the time 
            List<double> tempY = dt[1];//second is the RUs
                //two columns of course are identical in length
            List<double> tempXA = new List<double>();
            List<double> tempXD=new List<double>();
            this.C_Y = new List<double>();
            this.C_Y_Detach = new List<double>();
            for (int i=0; i < tempX.Count; i++)
            {
                if (tempX[i] < this.C_Duration_Attach)//attaching
                {
                    tempXA.Add(tempX[i]);
                    this.C_Y.Add(tempY[i]);
                }
                else //detaching
                {
                    tempXD.Add(tempX[i] - this.C_Duration_Attach);
                    this.C_Y_Detach.Add(tempY[i]);
                }
            }
            //C_Y = dt[2];
            this.C_X.Add(tempXA);
            this.C_X.Add(tempXD);
        }

        //member declaration
        List<double> C_Y_Detach;//for detach phase of data values
        double C_Conc;
        double C_Duration_Attach;//this is the cutoff value between attaching and detaching phases.
        
    }//end of class
}
