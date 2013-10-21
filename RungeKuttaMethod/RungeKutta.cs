using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace RungeKuttaMethod
{
    public delegate double FuctionDelegate(double t, double y);
    /// <summary>
    /// this is the fourth order RungeKutta Implementation
    /// </summary>
    public class RungeKutta
    {
        /// <summary>
        /// empty constructor
        /// </summary>
        public RungeKutta()
        {
        }
        /// <summary>
        /// the actual method implementing the RungeKutta method
        /// </summary>
        /// <param name="_fd">the derivate of the funtion to be numerical estimated</param>
        /// <param name="_input">the pair (t0, tmax) list, assume the initial values are the first one</param>
        /// <param name="_h">the step value</param>
        /// <param name="_initY">the intial value for the _input[1]</param>
        /// <returns>the y values for the (t,y)</returns>
        public static List<double> Solution(FuctionDelegate _fd, List<double> _input, double _h, double _initY)
        {
            List<double> ret = new List<double>();
            double k1,k2,k3,k4, currentY;

            ret.Add(_initY);
            for (int i = 1;i<=(_input[1]-_input[0])/_h ;i++ )
            {
                k1 = _fd(_input[0]+(i-1)*_h, ret[i-1]);
                k2 = _fd(_input[0] + (i-1) * _h + 0.5 * _h, ret[i - 1] + k1 * 0.5 * _h);
                k3 = _fd(_input[0] + (i-1) * _h + 0.5 * _h, ret[i - 1] + k2 * 0.5 * _h);
                k4 = _fd(_input[0] + (i-1) * _h + _h, ret[i - 1] + k3 * _h);

                currentY = ret[i - 1]+_h*(k1+2*k2+2*k3+k4)/6.0;

                ret.Add(currentY);
            }

            return ret;
        }
        /// <summary>
        /// the Runge-Kutta method for numerical solution
        /// </summary>
        /// <param name="_fd">the derivative function to be solved</param>
        /// <param name="_input">the independent variable, eg, x or t in many cases, this one lists all the points to be estimated,
        /// not just the bound values </param>
        /// <param name="_output">the ouput array, we assume this one has been initialized to have the equal length as the input, 
        /// with the initial value at the beginning, in many cases is zero. also this will be the output two </param>
        public static void Solution(FuctionDelegate _fd, List<double> _input, ref List<double> _output )
        {
            double k1, k2, k3, k4, currentY;
            //check whether the two arrays are the same.
            if (_input.Count != _output.Count)
                throw new System.Exception("the input and output is not set up correctly");

            for (int i = 1; i < _input.Count; i++)
            {
                double h = _input[i] - _input[i - 1];
                k1 = _fd(_input[i-1], _output[i - 1]);
                k2 = _fd(_input[i-1] + 0.5 * h, _output [i - 1] + k1 * 0.5 * h);
                k3 = _fd(_input[i-1] + 0.5 * h, _output[i - 1] + k2 * 0.5 * h);
                k4 = _fd(_input[i-1] + h, _output [i - 1] + k3 * h);

                currentY = _output[i - 1] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;

                _output[i]= currentY;
            }
        }

        /// <summary>DEPRICATED!!!!!!!!!!!!!
        /// this implementation is not right!!!!!, there is a H*fd()*1/2  in the term, just becareful.
        /// </summary>
        /// <param name="_fd"></param>
        /// <param name="_input"></param>
        /// <param name="_h"></param>
        /// <param name="_initY"></param>
        /// <returns></returns>
        public static List<double> SolutionH(FuctionDelegate _fd, List<double> _input, double _h, double _initY)
        {
            List<double> ret = new List<double>();
            double k1, k2, k3, k4, currentY;

            ret.Add(_initY);
            for (int i = 1; i <= (_input[1] - _input[0]) / _h; i++)
            {
                k1 = _fd(_input[0] + (i - 1) * _h, ret[i - 1]);
                k2 = _fd(_input[0] + (i - 1) * _h + 0.5 * _h, ret[i - 1] + k1 * 0.5 );
                k3 = _fd(_input[0] + (i - 1) * _h + 0.5 * _h, ret[i - 1] + k2 * 0.5 );
                k4 = _fd(_input[0] + (i - 1) * _h + _h, ret[i - 1] + k3 );

                currentY = ret[i - 1] + _h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;

                ret.Add(currentY);
            }

            return ret;
        }

    }//end of class
}//end of namespace.
