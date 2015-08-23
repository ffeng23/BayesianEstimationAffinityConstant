using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using Meta.Numerics;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;
using RungeKuttaMethod;

namespace BayesianEstimateLib
{
    /// <summary>
    /// this class is the one use Euler Method to run the numerical integration of SPR dynamics
    /// the reference for the dynamics is as follow,
    /// 1)Glaser, R. W. (1993). Antigen-antibody binding and mass transport by convection and diffusion to a surface: a two-dimensional computer model of binding and dissociation kinetics. Anal Biochem, 213(1), 152-161. doi: 10.1006/abio.1993.1399
    /// 2)Karlsson, R. (1994). Real-time competitive kinetic analysis of interactions between low-molecular-weight ligands in solution and surface-immobilized receptors. Anal Biochem, 221(1), 142-151. doi: 10.1006/abio.1994.1390
    /// In this one, we 
    /// </summary>
    public class NumericalIntegrationOfDynamics:SimulationSPR
    {
        /* don't allow empty constructor in order to keep the integrity of data
        /// <summary>
        /// empty constructer
        /// </summary>
        public NumericalIntegrationOfDynamics()
        {
            //empty constructer
        }
         * */
        /// <summary>
        /// the constructer with necessary parameters
        /// </summary>
        /// <param name="_ka">on rate constant, unit 1/m*1/s, 1E4 - 1E6</param>
        /// <param name="_kd">off rate constant, unit -/s, 1E-3 - 1E-6 (??)</param>
        /// <param name="conc">the concentration of analytes in the flow buffer, unit M (molar/l)</param>
        /// <param name="Rmax">the maximum Response Unit, unit arbitrary unit</param>
        /// <param name="duration">the association duration, unit second</param>
        /// <param name="kM">the diffusion rate, unit m/s, about 1E-5 in the BIAcore flow cell.</param>
        /// <param name="deltaT">time step to run numerical integration, the smaller the better, 1/100s or smaller</param>
        public NumericalIntegrationOfDynamics(double _ka, double _kd, double _conc, double _Rmax, double _r0,
                                                double _duration_attach, double _duration_detach, 
                                                double _kM=1E6, double _deltaT=0.01 ):base(_ka, _kd, _conc, _Rmax, _r0, 
                                                 _duration_attach, _duration_detach,  _kM, _deltaT)
        {
            //has been done in base.
        }
        /// <summary>
        /// same as the Euler scheme run_attach, but with RungeKutta method
        /// </summary>
        public override void run_Attach()
        {
            RungeKutta.Solution(this.DerivativeFunction_Attach, this._time_attach, ref this._ru_attach);
  
        }
        /// <summary>
        /// 
        /// same as the Euler scheme run_attach, but with RungKutta method
        /// </summary>

        public override void run_Detach()
        {
            _ru_detach[0] = this.SSPR_r0;
            RungeKutta.Solution(this.DerivativeFunction_Detach, this._time_detach, ref this._ru_detach);
        }

        /// <summary>Euler scheme
        /// run to integrate numerically to get the simulation
        /// the equation is as follow
        /// d[AB]/dt = kf*conc*([AB]max-[AB])-kr*[AB];
        /// kf = ka*kM/(kM+ka*([AB]max-[AB]));
        /// kr = kd*kM/(kM+ka*([AB]max-[AB]));
        /// d[AB]/dt=q(d[R]/dt)
        /// </summary>
        public void run_AttachEuler()
        {
            //_ru.Add(0);the _ru_attach has been initialized and added with all zeros in the base class.
            for( int i=0; ;i++)
            {
	        
	            double kf=_ka*_kM/(_kM+_ka*(_Rmax-_ru_attach[i]));
	            double kr=_kd*_kM/(_kM+_ka*(_Rmax-_ru_attach[i]));
	            double deltaR = kf*_conc*(_Rmax-_ru_attach[i])-kr*_ru_attach[i];
	            if(i>=_ru_attach.Count -1 )
	                {
                        break;
		                
	                }
                _ru_attach[i + 1] = deltaR * (_time_attach[i+1]-_time_attach[i]) + _ru_attach[i];
            }
        }
        /// <summary>Euler scheme
        /// still using the same equations as run_attach. see above, but with [conc]=0, starting at R0
        /// </summary>
        /// <param name="_R0"></param>
        public void run_DetachEuler()
        {
            
            _ru_detach[0]=this.SSPR_r0 ;
            
            //_ru.Add(0);
            for (int i = 0; ; i++)
            {

                double kf = _ka * _kM / (_kM + _ka * (_Rmax - _ru_detach[i]));
                double kr = _kd * _kM / (_kM + _ka * (_Rmax - _ru_detach[i]));
                double deltaR = 0 - kr * _ru_detach[i];
                if (i >= _ru_detach.Count - 1)
                {
                    break;

                }
                _ru_detach[i + 1] = deltaR * (_time_detach[i + 1] - _time_detach[i]) + _ru_detach[i];
            }
        }

        /// <summary>
        /// this is the function delegate for runge kutta method. this is the derivative funciton specifying the time differential relation
        /// </summary>
        /// <param name="t">time t</param>
        /// <param name="y">y value at time t, the RUs</param>
        /// <returns>the dy/dt value</returns>
        protected double DerivativeFunction_Attach(double _t, double _y)
        {
            double kf = _ka * _kM / (_kM + _ka * (_Rmax - _y));
            double kr = _kd * _kM / (_kM + _ka * (_Rmax - _y));
            double deltaR = kf * _conc * (_Rmax -_y) - kr * _y;

            return deltaR;
        }
        /// <summary>
        /// same as above, but for detaching phase
        /// </summary>
        /// <param name="_t"></param>
        /// <param name="_y"></param>
        /// <returns></returns>
        protected double DerivativeFunction_Detach(double _t, double _y)
        {
            double kf = _ka * _kM / (_kM + _ka * (_Rmax - _y));
            double kr = _kd * _kM / (_kM + _ka * (_Rmax - _y));
            double deltaR = 0 - kr * _y;

            return deltaR;
        }

        /// <summary>
        /// not implemented, don't use
        /// </summary>
        /// <param name="_params"></param>
        public override void setParameters(double[] _params)
        {
            throw new NotImplementedException();
        }


    }//end of class

}
