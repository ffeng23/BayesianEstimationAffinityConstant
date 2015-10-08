using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BayesianEstimateLib
{
    /// <summary>
    /// this is the Langumir model with analytical solution derived from SimulationSPR.
    /// </summary>
    public sealed class Langmuir:SimulationSPR 
    {
        public Langmuir(double _ka, double _kd, double _conc, double _Rmax, double _r0,
                                                double _duration_attach, double _duration_detach, double _deltaT=0.01 )
            :base(_ka,  _kd,  _conc, _Rmax, _r0, _duration_attach,  _duration_detach, 0, _deltaT)
            {
                //doing nothing
            }
        public Langmuir():base()
        {
            //everything has been initialized to -1 or null
        }

        /// <summary>
        /// Langmuir model
        /// Rt=ka*Conc*Rmax/(kd+ka*Conc)*(1-exp(-(kd+ka*Conc)t)
        /// </summary>
        public override void run_Attach()
        {
            //_ru.Add(0);
            for (int i = 0; ; i++)
            {
                if(i >= _ru_attach.Count )
                {
                    break;

                }
                _ru_attach[i] = _ka *_conc *_Rmax /(_kd +_ka *_conc)*(1-Math.Exp(-1*(_kd +_ka *_conc)*_time_attach[i]));
            
                
                //_ru_attach[i + 1] = deltaR * (_time_attach[i + 1] - _time_attach[i]) + _ru_attach[i];
            }
        }
        /// <summary>
        /// langmuir model, starting at R0
        /// </summary>
        /// <param name="_R0"></param>
        public override void run_Detach()
        {
            if (this.SSPR_r0 > 0)
                _ru_detach[0] = this.SSPR_r0;
            else
                _ru_detach[0] = this._ru_attach[_ru_attach.Count() - 1];
            //_ru.Add(0);
            for (int i = 0; ; i++)
            {
                if (i >= _ru_detach.Count )
                    {
                        break;

                    }
                _ru_detach[i] = this._ru_detach[0]  * Math.Exp(-1*_time_detach[i] *_kd) ;
                
                
            }
        }
        /// <summary>
        /// set the parameter according to the input _params array
        /// 0. ka; 1,kd; 2-conc, 3-Rmax, 4-R0_AB,
        ///      5-association duration, 6-association duration, 7-deltaT
        ///     So far DON"T have kM in there now.
        ///     here we also allow variable length of params, but also
        ///     allow default value of deltaT, duration of association or detach.
        /// </summary>
        /// <param name="_params">0. ka; 1,kd; 2-conc, 3-Rmax, 4-R0_AB,
        ///      5-association duration, 6-association duration, 7-deltaT</param>
        public override void setParameters(double[] _params)
        {
            //here we allow variable length
            if (_params.Count() < 5) //5 parameters are the minimum for repeated calling of established object
            //in this case, we don't have to repeat the filling time arrays etc.
            //probability 4 is the minimum, but we do this 5 for now.
            {
                throw new Exception("the input parameter array is not valid. doesn't contain enough elements");
            }
            this._ka = _params[0];
            this._kd = _params[1];
            
            this._conc = _params[2];
            this._Rmax = _params[3];
            this.SSPR_r0 = _params[4];
            //this._R0_AB_Star = _params[7];

            bool updateTimeArrays = false;
            if (_params.Count() >= 6 && _params[5] > 0)
            {
                updateTimeArrays = true;
                this._duration_attach = _params[5];
            }
            if (_params.Count() >= 7 && _params[6] > 0)
            {
                updateTimeArrays = true;
                this._duration_detach = _params[6];
            }
            if (_params.Count() >= 8 && _params[7] > 0)
            {
                updateTimeArrays = true;
                this._deltaT = _params[7];
            }

            if (updateTimeArrays)
            {
                _fillTimeArrays();
            }
        }
    }//end of class
}
