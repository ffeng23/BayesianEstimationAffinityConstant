using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace RungeKuttaMethod
{
    public delegate double FuctionDelegate(double t, double y);
    public delegate List<double> FunctionDelegates(double t, List<double> ys);
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
        /// <summary>
        /// the RK 4th order to do numerical integration of functions containing MULTIPLE variables
        /// Y={y1,y2...yi}
        /// Y=f(Y,t)
        /// it is important to know that t is a scalar and Y is a vector for each steps of small dt
        /// 
        /// </summary>
        /// <param name="_fds">function delegate takes t and yi at time t</param>
        /// <param name="_input">input is the a list of t with small step dt</param>
        /// <param name="_output">this is both the input and output. it must have first element define and based on
        /// this initial Yi0, we will generate Yij at each input t. ALSO very importantly output array 
        /// has to be initialize with the equal number of t.
        /// </param>
        public static void Solution(FunctionDelegates _fds, List<double > _input, ref List<List <double> >_output)
        {
            List<double> k1; List<double> k2; List<double> k3;List<double> k4;
            //List<double> currentYs;
            double h;
            //check whether the two arrays are the same.
            if (_input.Count != _output.Count)
                throw new System.Exception("the input and output is not set up correctly");

            List<double> temp_array=new List<double>();

            for (int i = 1; i < _input.Count; i++)
            {
                h = _input[i] - _input[i - 1];
                for (int j = 0; j < _output[0].Count; j++)//
                {
                    temp_array=_output[i-1];
                    k1 = _fds(_input[i - 1], temp_array);
                    temp_array = _UpdateArray(_output[i - 1], k1, 0.5 * h);
                    k2 = _fds(_input[i - 1] + 0.5 * h, temp_array);
                    temp_array = _UpdateArray(_output[i - 1], k2, 0.5 * h);
                    k3 = _fds(_input[i - 1] + 0.5 * h, temp_array);
                    temp_array = _UpdateArray(_output[i - 1], k3,  h);
                    k4 = _fds(_input[i - 1] + h, temp_array);

                    temp_array = _UpdateArray_4RK(_output[i - 1], k1 , k2 ,k3,k4, h);

                    _output[i]= temp_array;
                }
            }
        }

        /// <summary>
        /// array operation, used to calculate _output[i-1]+k1*0.5*h)
        /// </summary>
        /// <param name="_a"></param>
        /// <param name="_k"></param>
        /// <param name="_fact"></param>
        /// <returns></returns>
        protected static List<double> _UpdateArray(List<double> _a, List<double> _k, double _fact)
        {
            List<double> ret_a = new List<double>();

            for (int i = 0; i < _a.Count; i++)
            {
                ret_a.Add(_a[i] + _k[i] * _fact);
            }
            return ret_a;
        }
        /// <summary>
        /// do the calculation as currentY = _output[i - 1] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
        /// </summary>
        /// <param name="_a"></param>
        /// <param name="_k1"></param>
        /// <param name="_k2"></param>
        /// <param name="_k3"></param>
        /// <param name="_k4"></param>
        /// <returns></returns>
        protected static List<double> _UpdateArray_4RK(List<double> _a, List<double> _k1,
            List<double> _k2, List<double> _k3, List<double> _k4, double _h)
        {
            //currentY = _output[i - 1] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
            List<double> currentY = new List<double>();
            for (int i = 0; i < _a.Count; i++)
            {
                currentY.Add(_a[i]+_h*(_k1[i]+2*_k2[i]+2*_k3[i]+_k4[i])/6.0);
            }
            return currentY;
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
