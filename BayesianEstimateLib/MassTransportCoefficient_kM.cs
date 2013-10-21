using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BayesianEstimateLib
{
    public class MassTransportCoefficient_kM
    {
        //**************constant variable used for calculate diffusion coefficeint. we will need input of temperature and MW.
        static double k = 1.381E-23; //boltzmann's constant, unit J/K
        //static const double tc; //temperature in centigrade. will be inputted 
        static double TConst = 273.15;//the absolute temperature constant;
        static double v = 7.0E-4; //specific density, for protein, normally is 7.0-7.55E-4, unit m^3/kg
        static double Na = 6.022E+23; //Avogardro's number, unit /mol
        static double eta = 0.001; //viscosity, unit kg/(s*m)
        static  double flf0 = 1.2; //friction factor 1.2 a reasonable guess for most globular protein
        //static const double MW; //molecular weight, an input, unit kg/mol (kDa) for IgG

        /// <summary>
        /// 
        /// </summary>
        /// <param name="_tc">degree</param>
        /// <param name="_mw">kg/mol or kDa</param>
        /// <returns>m^2/s</returns>
        public static double CalculateDiffusionCoefficient(double _tc, double _mw)
        {
            double D=0;
            D=k*(_tc+TConst)/( 6*Math.PI*Math.Pow(3*_mw*v/(4*Math.PI *Na),1.0/3.0)*eta*flf0);

            return D; 
        }

        // for kM
        // double L1, L2, h, L, w; //all are input as the SPR dection area
        //double F; //flow rate as input
        static  double G = 1000; // this is magic number to convert between RU and attached Mass (g), unit RU* mm^2/ng

        /// <summary>
        /// this one is used to calculate the kM in unit of m/s
        /// </summary>
        /// <param name="_tc">degrade</param>
        /// <param name="_mw">kDa</param>
        /// <param name="_F">unit ul/min</param>
        /// <param name="_L1">unit mm</param>
        /// <param name="_L2">unit mm</param>
        /// <param name="_h">unit mm</param>
        /// <param name="_w">unit mm</param>
        /// <returns>m/s</returns>
        public static double CalculatekM(double _tc, double _mw, double _F, double _L1, double _L2, double _h, double _w, double _D)
        {
            double kM = 0;
            double CkM = 1.47 * (1 - Math.Pow((_L1 / _L2), 2.0 / 3.0)) / (1 - _L1 / _L2);
            //double D = CalculateDiffusionCoefficient(_tc, _mw);
            //now we need to convert to standard units
            double _L1_prime = _L1 / 1000;
            double _L2_prime = _L2 / 1000;
            double _h_prime = _h / 1000;
            double _w_prime = _w / 1000;
            double _F_prime = _F / 1E9 / 60;
            kM = CkM * Math.Pow(_D * _D * _F_prime / (_h_prime * _h_prime * _w_prime * _L2_prime), 1.0 / 3.0);//now this is in unit of m/s
            
            return kM;
        }
        /// <summary>
        /// 
        /// </summary>
        ///<param name="_mw">g/mol</param>
        ///<param name="_kM">m/s</param>
        /// <returns>RU/(M*s)</returns>
        public static double CalculatekM_prime( double _mw, double _kM)
        {
            double kM_prime = 0;
            //double kM = CalculatekM(_tc, _mw, _F, _L1, _L2, _h, _w);
            double G_prime = G * Math.Pow((1E-2), 2) / 1E-9; //now it is in RU*(dm)^2/g 
            kM_prime = _kM * G_prime * (_mw * 1000) *10; //10 is used to conver m to dm.
            return kM_prime; //now it is RU*dm^3/mol/s or RU/(M*s)
        }


    }//end of 
}
