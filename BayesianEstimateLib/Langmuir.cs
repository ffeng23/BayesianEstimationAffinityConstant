using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BayesianEstimateLib
{
    public sealed class Langmuir:SimulationSPR 
    {
        public Langmuir(double _ka, double _kd, double _conc, double _Rmax, double _r0,
                                                double _duration_attach, double _duration_detach, double _deltaT=0.01 )
            :base(_ka,  _kd,  _conc, _Rmax, _r0, _duration_attach,  _duration_detach, 0, _deltaT)
            {
                //doing nothing
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

            _ru_detach[0] = this.SSPR_r0 ;

            //_ru.Add(0);
            for (int i = 0; ; i++)
            {
                if (i >= _ru_detach.Count )
                    {
                        break;

                    }
                _ru_detach[i] = this.SSPR_r0  * Math.Exp(-1*_time_detach[i] *_kd) ;
                
                
            }
        }

        public override void setParameters(double[] _params)
        {
            throw new NotImplementedException("I am not done yet");
        }
    }//end of class
}
