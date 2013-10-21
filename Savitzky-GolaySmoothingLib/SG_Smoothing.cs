using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
//using Savitzky_GolaySmoothingLib;

/*************************README***************
 *This is the well known method to smooth the experimental data by convalution. 
 * based on 3 References
 * 1.	Savitzky, A. and M.J.E. Golay, Smoothing and Differentiation of Data by Simplified Least Squares Procedures. Analytical Chemistry, 1964. 36(8): p. 1627-1639.
 * 2.	Gorry, P.A., General least-squares smoothing and differentiation by the convolution (Savitzky-Golay) method. Analytical Chemistry, 1990. 62(6): p. 570-573.
   3.	Madden, H.H., Comments on the Savitzky-Golay convolution method for least-squares-fit smoothing and differentiation of digital data. Analytical Chemistry, 1978. 50(9): p. 1383-1386. 
 *
 * Basically, we are running the polynomial regression through the datapoints and get the estimate smooth curve and derivatives.
 * here we implement the general method as in Ref 2. 
 */

namespace Savitzky_GolaySmoothingLib
{
    /// <summary>
    /// 
    /// </summary>
    public class SG_Smoothing
    {
        //empty constructor
        //this is enabled
        /*private SG_Smoothing()
        {
            //empty constructor
        }*/

        /// <summary>
        /// constructor
        /// </summary>
        /// <param name="_frameSize">window size is normally 9 (4*2+1)</param>
        /// <param name="_degreeOfPolynomial">polynomial order of 4 by default</param>
        public SG_Smoothing(int _frameSize=4, int _degreeOfPolynomial=4)
        {
            //setting the data
            C_degreeOfPolynomial = _degreeOfPolynomial ;
            C_frameSize = _frameSize;
        }

        public void SetParameters(int _frameSize , int _degreeOfPolynomial )
        {
            //setting the data
            C_degreeOfPolynomial = _degreeOfPolynomial;
            C_frameSize = _frameSize;
        }
        /// <summary>
        /// this is the working function to generate the Coeffecients matrix/dictionary, keyed by the position of points
        /// the size(count) of dictionary is determined by the frame size (frameSize), 2m+1. the posittionn is the "t" in tthe
        /// paper, (-m<= t <=m)
        /// </summary>
        /// <param name="OrderOfDerivative">order of derivatives</param>
        /// <returns></returns>
        public Dictionary<int, List<List<BigInteger >>> GenerateCoefficients(int _orderOfDerivative)
        {
            Dictionary<int, List<List<BigInteger >>> cfm = new Dictionary<int, List<List<BigInteger >>>();
            for (int t = -1 * this.C_frameSize; t <= this.C_frameSize; t++) //for each position
            {
                //we need to get all the index coefficients
                List<List<BigInteger >> cfm_t = new List<List<BigInteger >>();
                for (int i = -1 * this.C_frameSize; i <= this.C_frameSize; i++)
                {
                    //get the coefficient
                    cfm_t.Add( Weight(i, t, C_frameSize, C_degreeOfPolynomial, _orderOfDerivative));
                   
                } 
                cfm.Add(t, cfm_t);
            }

            return cfm;
        }

        public List<Double> Smooth(int _orderOfDerivative, List<Double> _data)
        {
            List<Double> output = new List<Double>();

            //first get coefficient matrix
            Dictionary<int,List<List<BigInteger >>> cfm=GenerateCoefficients(_orderOfDerivative );
            List<List<BigInteger >> cfm_t;
            
            if (_data.Count < 2 * C_frameSize)
            {
                Console.WriteLine("Too few data points");
                throw new System.Exception("too few data points");
                    
            }
            for (int i = 0; i < _data.Count; i++)
            {
                List<Double> dataPoints = new List<Double>(2 * C_frameSize + 1);
                int startingIndex;
                //taking care the initial points first
                if (i < C_frameSize)
                {
                    //get initial points coeffiecient matrix
                    cfm_t = cfm[i - C_frameSize];
                    startingIndex = 0;

                }
                else
                {
                    if (_data.Count - i - 1 < C_frameSize)
                    {//the last initial points
                        cfm_t = cfm[C_frameSize - (_data.Count - i - 1)];
                        startingIndex = _data.Count - 2*C_frameSize-1;
                    }
                    else //the middle ones using t=0 /poisition is centered
                    {
                        cfm_t = cfm[0];
                        startingIndex = i - C_frameSize;
                    }
                }//else outer
                //get datapoint array
                output.Add(0);
                for (int j = -1 * C_frameSize; j <= C_frameSize; j++)
                {
                    double tem=BigIntegerDivsion (cfm_t[j+C_frameSize ][0],cfm_t[j+C_frameSize][1]);
                    output[i] += _data[startingIndex + j + C_frameSize] * tem ;
                }

            }//end of for loop outer

            return output;
        }
        public static double BigIntegerDivsion(BigInteger bi1, BigInteger bi2)
        {
            int finalSign = 1;
            if ((bi1.Sign == -1 && bi2.Sign == 1) || (bi1.Sign == 1 && bi2.Sign == -1))
            {
                finalSign = -1;
            }
            double te = BigInteger.Log10(BigInteger.Abs(bi1)) - BigInteger.Log10(BigInteger.Abs(bi2));
            double teE=Math.Pow(10,te);
            return finalSign * teE;
        }

        /// <summary>
        /// calculate the grand polynomial weights at i. framesize and degree of polynomial has been set through the class
        /// </summary>
        /// <param name="_index">i </param>
        /// <param name="_OrderOfDerivative">s</param>
        /// <returns>return the numerator and denominator seperately as an array</returns>
        public static List<BigInteger> GrandPolynomial(int _index, int _orderOfDerivative, int _degreeOfPolynomial, int _frameSize)
        {
            List<BigInteger > retList = new List<BigInteger>(2);//return array, contains first nominator and second denominator
            //retList.Add(0); retList.Add(0);

            BigInteger numerator, denominator;

            if (_degreeOfPolynomial > 0)
            {
                List<BigInteger > GP_Kminus1 = GrandPolynomial(_index, _orderOfDerivative, _degreeOfPolynomial - 1, _frameSize);
                List<BigInteger > GP_Kminus1_Sminus1 = GrandPolynomial(_index, _orderOfDerivative - 1, _degreeOfPolynomial - 1, _frameSize);
                List<BigInteger> GP_Kminus2 = GrandPolynomial(_index, _orderOfDerivative, _degreeOfPolynomial - 2, _frameSize);

                //numerator
                BigInteger  numerator1 = (4 * _degreeOfPolynomial - 2) * (_index * GP_Kminus1[0] * GP_Kminus1_Sminus1[1] + _orderOfDerivative * GP_Kminus1_Sminus1[0] * GP_Kminus1[1]);
                numerator1 = numerator1 * GP_Kminus2[1];
                BigInteger numerator2 = (_degreeOfPolynomial - 1) * (2 * _frameSize + _degreeOfPolynomial) * GP_Kminus2[0];
                numerator2 = numerator2 * GP_Kminus1[1] * GP_Kminus1_Sminus1[1];

                numerator = numerator1 - numerator2;

                //denominator
                denominator = _degreeOfPolynomial * (2 * _frameSize - _degreeOfPolynomial + 1) * GP_Kminus1[1] * GP_Kminus1_Sminus1[1] * GP_Kminus2[1];

                retList.Add(numerator);
                retList.Add(denominator);
                SimplifyNumDenom(retList);
            }
            else
            {
                if (_degreeOfPolynomial == 0 && _orderOfDerivative == 0)
                {
                    retList.Add(1); retList.Add(1);
                }
                else
                {
                    retList.Add(0); retList.Add(1);
                }
            }

            //SimplifyNumDenom(retList);
            return retList;
        }

        public static BigInteger GeneralFactorial(int a, int b)
        {
            BigInteger  gf = 1;
            for (int j = a - b + 1; j <= a; j++)
            {
                gf *= j;
            }

            return gf;
        }

        public static List<BigInteger> Weight(int _index, int _positionOfpoints, int _frameSize, int _degreeOfPolynomial, int _orderOfderivative)
        {
            List<BigInteger> retList = new List<BigInteger >(2);
            BigInteger numerator=0, denominator=1;

            for (int k = 0; k <= _degreeOfPolynomial; k++)
            {
                List<BigInteger> GP_index=GrandPolynomial(_index, 0, k, _frameSize);
                List<BigInteger> GP_position=GrandPolynomial(_positionOfpoints , _orderOfderivative,k,_frameSize);

                numerator = numerator * GeneralFactorial(2 * _frameSize + k + 1, k + 1) * GP_index[1] * GP_position[1];
                numerator = numerator + denominator * (2 * k + 1) * GeneralFactorial(2 * _frameSize, k) * GP_index[0] * GP_position[0];

                
                denominator = denominator * GP_index[1] * GP_position[1] * GeneralFactorial(2 * _frameSize + k + 1, k + 1);
                SimplifyNumDenom(ref numerator, ref denominator);
                //Console.WriteLine("after k=" + k + ", numerator is " + numerator + "; deno is " + denominator);
            }

            retList.Add(numerator); retList.Add(denominator);
            //SimplifyNumDenom(retList);
            return retList;
        }

        /// <summary>
        /// used to simply the coFactor of numerator and denominator, it works in situ.
        /// </summary>
        /// <param name="x">intput list as [0]numerator and [1]denominator</param>
        public static void SimplifyNumDenom(List<BigInteger> x)
        {


            if (x[0] == 0)
            {
                x[1] = 1;
            }
            List<BigInteger> primeNumbers = new List<BigInteger> { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997 };
            for (int i = 0; i < primeNumbers.Count; i++)
            {
                if (x[0] < primeNumbers[i]&&x[0]>-1*primeNumbers[i])
                {
                    break;
                }
                while (x[0] / primeNumbers[i] * primeNumbers[i] == x[0] && x[1] / primeNumbers[i] * primeNumbers[i] == x[1])
                {
                    x[0] /= primeNumbers[i];
                    x[1] /= primeNumbers[i];
                }
            }
        }
        /// <summary>
        /// a version that takes in reference of two 
        /// </summary>
        /// <param name="_n"></param>
        /// <param name="_d"></param>
        public static void SimplifyNumDenom(ref BigInteger _n, ref BigInteger _d)
        {


            if (_n == 0)
            {
                _d = 1;
            }
            List<int> primeNumbers = new List<int> { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997 };
            for (int i = 0; i < primeNumbers.Count; i++)
            {
                if (_n < primeNumbers[i]&&_n>-1*primeNumbers[i])
                {
                    break;
                }
                while (_n / primeNumbers[i] * primeNumbers[i] == _n && _d / primeNumbers[i] * primeNumbers[i] == _d)
                {
                    _n /= primeNumbers[i];
                    _d /= primeNumbers[i];
                }
            }
        }

        //****************member declaration
        private int C_degreeOfPolynomial; //k
        //private int C_degreeOfDerivative; //s
        //private int C_datapointIndex; // i in the ref paper
        //private int C_positionOfPolynomial; //t 
        private int C_frameSize; //m in the paper, total data points is 2m+1;


    }//end of class

     
}
