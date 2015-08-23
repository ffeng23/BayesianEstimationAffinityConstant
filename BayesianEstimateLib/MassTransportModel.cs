using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BayesianEstimateLib
{
    /// <summary>
    /// this class is the same as the numericalIntegrationOfDynamics. use that one please 
    /// </summary>
    public sealed class MassTransportModel:SimulationSPR 
    {
        public MassTransportModel(double _ka, double _kd, double _conc, double _Rmax, double _r0,
                                                double _duration_attach, double _duration_detach, double _kM = 1E6, double _deltaT = 0.01)
            :base(_ka, _kd, _conc, _Rmax, _r0,_duration_attach,_duration_detach,_kM, _deltaT)
        {
            //empty
        }

        public override void run_Attach()
        {
            //_ru.Add(0);
            for (int i = 0; ; i++)
            {

                
                double deltaR = _kM * ( _conc*_ka*(_Rmax-_ru_attach[i])  - _kd * _ru_attach[i])/(_kM+_ka*(_Rmax-_ru_attach[i]));
                if (i >= _ru_attach.Count - 1)
                {
                    break;

                }
                _ru_attach[i + 1] = deltaR * (_time_attach[i + 1] - _time_attach[i]) + _ru_attach[i];
            }
        }
        /// <summary>
        /// still using the same equations as run_attach. see above, but with [conc]=0, starting at R0
        /// </summary>
        /// <param name="_R0"></param>
        public override void run_Detach()
        {

            _ru_detach[0] = this.SSPR_r0;

            //_ru.Add(0);
            for (int i = 0; ; i++)
            {

                double deltaR = _kM * (0 - _kd / _ka * _ru_detach[i] / (_Rmax - _ru_detach[i])); ;
                if (i >= _ru_detach.Count - 1)
                {
                    break;

                }
                _ru_detach[i + 1] = deltaR * (_time_detach[i + 1] - _time_detach[i]) + _ru_detach[i];
            }
        }
        public override void setParameters(double[] _params)
        {
            //don't use this one.
            throw new Exception("don't use me yet, since I am not done");
        }
    }//end of class
}
