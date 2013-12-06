using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace NelderMeadMethod
{
    public delegate double LogDistributionFuctionDelegate(double x, double functionNormConstant);
    
    /// <summary>
    /// this is the nelderMead algorithm in one dimension case. so far this is used to find the nozero value of any given function
    /// this will use the log value of the function and going up to find the nonzero value of the function. 
    /// </summary>
    public class NelderMead1D
    {   
        public  const  double EPSILON=1E-100;
        public const double LOG_LIMIT = -740;
        public NelderMead1D()
        {
        }

        /// <summary>
        /// this is the working algorithm to find the nonzero value applying N-M algorithm
        /// </summary>
        /// <param name="_logfunc">the log func to evaluated</param>
        /// <param name="_x0">the initial value for the algorithm</param>
        /// <param name="_x1">the initial value for the algorithm</param>
        /// <param name="_lower">the _lower bound for the domain</param>
        /// <param name="_upper">the _upper bound for the domain</param>
        /// <param name="_logLimit">the limit in log format that will need to reach to be nonzero value, so this means
        /// if we reach this value (of fucntion evaluation), we are basically done with the job. The reason we need this is 
        /// that we are using log function, and log(0) is not allowed, so we need to specify what is the value that is close to
        /// zero. 
        /// </param>
        /// <param name="_functionNormConstant">we need this one for some of the function that to normalize the value, in case 
        /// this _LogLimit is be not good. Log(x) is define by the C# to be -700., but sometimes, this is not very close enough to
        /// Zero. So we need to add this normal constant to scale it and we still are able to run function.
        /// </param>
        /// <returns></returns>
        public static double FindNonZeroValue(LogDistributionFuctionDelegate  _logfunc, double _x0, double _x1, 
            double _lower, double _upper, double _logLimit, ref double _functionNormConstant )
        {
            double xb;
            while (true)
            {
                double f0 = _logfunc(_x0, _functionNormConstant);
                if (f0 > _logLimit)
                {
                    return _x0;
                }
                double f1 = _logfunc(_x1, _functionNormConstant);
                if (f1 > _logLimit)
                {
                    return _x1;
                }
                //Console.WriteLine("\t\tIn foundNonZeroValue:x0 is " + _x0 + "; x1 is " + _x1 + ";fo is " + f0 + ";f1 is " + f1);
                xb = FindNonZeroValueKernelLarger(_logfunc, _x0, _x1, _lower, _upper, _logLimit, f0, f1, _functionNormConstant, C_MaximumIterations);
                double fxb = _logfunc(xb, _functionNormConstant);
                if (fxb<=_logLimit)
                {
                    _functionNormConstant = 300.0 - fxb;
                    
                }
                else
                    return xb;
            }
                
            
        }//end of function fidnNonzerovalue

        private static double FindNonZeroValueKernelLarger(LogDistributionFuctionDelegate  _logfunc, double _x0, double _x1, double _lower, double _upper,
            double _logLimit, double _f0, double _f1, double _functionNormConstant, int _maxIters)
        {
            double f0=_f0;
            
            double f1=_f1;
            double epsilon;

            //else if we are here, we need to run the N-M algorithm

            double fb, fw, xb, xw;//f best, f worst, x best, x worst
            double xCentroid;
            double xn, fn;//the new x, also used as fr, xr,the reflection point
            fb=f0;//so far, neither f0 nor f1 is good enough to quit, we start with either one to start the while loop
            xb = _x0;
            fn = f1;
            xn = _x1;

            double count = 0;
            while(true)  //keep going
            {
                count++;
                //Console.WriteLine("looping steps:" + count);
                double xtemp, ftemp;
                fw = fn; xw = xn;
                //ordering
                if (fn > fb) //swapping
                {
                    xtemp = xb; ftemp = fb;
                    fb = fn; xb = xn;
                    fw = ftemp; xw = xtemp;
                }
                epsilon = xb - xw;
                if ((epsilon < 1E-20 && epsilon > -1E-20) || (xb != 0 && (epsilon / xb < 1E-6 && epsilon / xb > -1E-6)))
                {
                    return xb;
                }
                if(fb>_logLimit) //we are done
                    break;
                if (_maxIters  != -1 && count > _maxIters)
                {
                    return xb;
                    //break;
                }
                /*else
                {
                    fb = f1; fw = f0;
                    xb = _x1; xw = _x0;
                }*/

                //here we need to consider the special boundary conditions, first
                //if xb is the boundary(_lower or upper), then we simply shrink the points towards the boundary
                if (
                    ((xb - _lower) < EPSILON && (xb - _lower)  > -1 * EPSILON) ||
                    ((xb - _upper)  < EPSILON && (xb - _upper)  > -1 * EPSILON)
                   )
                {
                    //we can simply try one more step to testing whether the reflection of this best point over the worst point will getting better???
                    xCentroid = xw;
                    xn = xCentroid + C_alpha * (xCentroid-xb); //reverse reflection over xw;
                    fn=_logfunc(xn, _functionNormConstant );
                    if (fn > fb)  //the reverse reflection gives a better point over the best
                    {
                        //here we simply replace the  points
                        xb=xn;
                        fb = fn;
                        xn = xw;
                        fn = fw;
                    }
                    else //the reverse reflection in above didnot becme better, so we shrink here.
                    {
                        xn = xb + C_sigma * (xw - xb);
                        fn = _logfunc(xn ,_functionNormConstant );
                    }
                    continue;
                }

                //if we are here, it means the boundaries (_lower, _upper) are not the best ones, so we can try to do the reflection now
                //doing reflection, xn=xcentroid. 
                //it is possible that the xw could be boundaries. but the reflection is going away from xw/boundaries in this case, so we are fine.
                xCentroid = xb;
                xn = xCentroid + C_alpha * (xCentroid - xw);

                //now we still have to test, whether this new ones are out of boundaries
                //it is possible that the either of xb, xw is not boundaries, but xb is closer to boundries, so xn ends up with being out of boundaries, 
                //so we have to check.
                if (xn <= _lower)
                {
                    xn = _lower;
                }
                if (xn >= _upper)
                {
                    xn = _upper;
                }
                fn = _logfunc(xn, _functionNormConstant);

                //we should never accept fn in here in this 1-d case, since fn is never between fs and fb, since fs and fb in this case is the same.
                

                //check what to do next
                if (fn < fb)  //so this mean we need to contract, in this case fn<fsecond for other cases, now we need to contract
                { //contract inside or outside???
                    double xc, fc;
                    if (fn > fw) //contract outside. this is better
                    {
                        xc = xCentroid + C_beta * (xn - xCentroid);
                        fc = _logfunc(xc, _functionNormConstant);
                        if (fc > fn)
                        {
                            fn = fc;
                            xn = xc;
                            continue;
                        }
                        //esle will go to do shrink
                    }
                    else //contract inside????
                    {
                        xc = xCentroid + C_beta * (xw - xCentroid);
                        fc = _logfunc(xc, _functionNormConstant);
                        if (fc > fn)
                        {
                            fn = fc;
                            xn = xc;
                            continue;
                        }
                        //else will go to do shrink

                    }
                }
                else //this means fn>fb, so we are doing expansion
                {
                    //compute the expansion point now,
                    double xe, fe;
                    xe = xCentroid + C_gamma * (xn - xCentroid);
                    fe = _logfunc(xe, _functionNormConstant);
                    if (fe >= fn)//will do the expansion
                    {
                        //but need to check whether we are out boundaries,
                        if (xe < _lower)
                        {
                            xe = _lower;
                            fe = _logfunc(xe, _functionNormConstant);
                            if (fe >= fn)//do the expansion
                            {
                                fn = fe;
                                xn = xe;
                                continue;
                            }
                            else //accept the fn xn as the reflection point,
                            {
                                continue;
                            }
                        }
                        if (xe > _upper)
                        {
                            xe = _upper;
                            fe = _logfunc(xe,_functionNormConstant);
                            if (fe >= fn)//do the expansion
                            {
                                fn = fe;
                                xn = xe;
                                continue;
                            }
                            else //accept the fn xn as the reflectio points
                                continue;
                        }
                        fn = fe;
                        xn = xe;
                        continue;
                    }
                    else //fn=fn, fn=xr; accept the reflection point
                    {
                        continue;
                    }
                }

                //if we are here, it means we need to run shrinkage.
                xn = xb + C_sigma * (xw - xb);
                fn = _logfunc(xn,_functionNormConstant);

            }//end of outer while

            return xb;
        }

        /// <summary>
        /// here this one is calling the kernels, then we use this one search the value bigger than the _logLimit in order to define the Envelop function tails at the lower end
        /// </summary>
        /// <param name="_logfunc"></param>
        /// <param name="_x0"></param>
        /// <param name="_x1"></param>
        /// <param name="_lower"></param>
        /// <param name="_upper"></param>
        /// <param name="_logLimit">the value that the function to search for a bigger one than this</param>
        /// <returns></returns>
        public static double FindBigValue(LogDistributionFuctionDelegate _logfunc, double _x0, double _x1, double _lower, double _upper
            /*,double _logLimit*/, double _functionNormConstant) //the _logLimit in this case is the f0=_logFunc(x0), always using f0/x0
        {
            double f0 = _logfunc(_x0, _functionNormConstant );

            double f1 = _logfunc(_x1, _functionNormConstant);

            double xb = FindNonZeroValueKernelLarger (_logfunc, _x0, _x1, _lower, _upper, f0, f0, f1,_functionNormConstant, -1/*maxiters*/);
            return xb;
        }

        /*
        public static double FindSmallValue(LogDistributionFuctionDelegate _logfunc, double _x0, double _x1, double _lower, double _upper
            , double _functionNormConstant)//the _logLimit in this case is the f0=_logFunc(x0)
        {
            double f0 = _logfunc(_x0, _functionNormConstant);

            double f1 = _logfunc(_x1, _functionNormConstant);

            double xb = FindNonZeroValueKernelSmaller(_logfunc, _x0, _x1, _lower, _upper, f0, f0, f1, _functionNormConstant,-1);
            return xb;
        }*/

        /*private static double FindNonZeroValueKernelSmaller(LogDistributionFuctionDelegate _logfunc, double _x0, double _x1, double _lower, double _upper,
            double _logLimit, double _f0, double _f1)
        {
            double f0 = _f0;
            double f1 = _f1;

            double fb, fw, xb, xw;//f best, f worst, x best, x worst
            double xCentroid;
            double xn, fn;//the new x, also used as fr, xr,the reflection point
            fb = f0;//so far, neither f0 nor f1 is good enough to quit, we start with either one to start the while loop
            xb = _x0;
            fn = f1;
            xn = _x1;

            double count = 0;
            while (true)  //keep going
            {
                count++;
                Console.WriteLine("looping steps:" + count);
                double xtemp, ftemp;
                fw = fn; xw = xn;
                //ordering
                if (fn < fb) //swapping
                {
                    xtemp = xb; ftemp = fb;
                    fb = fn; xb = xn;
                    fw = ftemp; xw = xtemp;
                }
                if (fb < _logLimit) //we are done
                    break;
         
                //else
                //{
                //    fb = f1; fw = f0;
                //    xb = _x1; xw = _x0;
                //}

                //here we need to consider the special boundary conditions, first
                //if xb is the boundary(_lower or upper), then we simply shrink the points towards the boundary
                if (
                    ((xb - _lower) < EPSILON && (xb - _lower) > -1 * EPSILON) ||
                    ((xb - _upper) < EPSILON && (xb - _upper) > -1 * EPSILON)
                   )
                {
                    //we can simply try one more step to testing whether the reflection of this best point over the worst point will getting better???
                    xCentroid = xw;
                    xn = xCentroid + C_alpha * (xCentroid - xb); //reverse reflection over xw;
                    fn = _logfunc(xn);
                    if (fn < fb)  //the reverse reflection gives a better point over the best
                    {
                        //here we simply replace the  points
                        xb = xn;
                        fb = fn;
                        xn = xw;
                        fn = fw;
                    }
                    else //the reverse reflection in above didnot becme better, so we shrink here.
                    {
                        xn = xb + C_sigma * (xw - xb);
                        fn = _logfunc(xn);
                    }
                    continue;
                }

                //if we are here, it means the boundaries (_lower, _upper) are not the best ones, so we can try to do the reflection now
                //doing reflection, xn=xcentroid. 
                //it is possible that the xw could be boundaries. but the reflection is going away from xw/boundaries in this case, so we are fine.
                xCentroid = xb;
                xn = xCentroid + C_alpha * (xCentroid - xw);

                //now we still have to test, whether this new ones are out of boundaries
                //it is possible that the either of xb, xw is not boundaries, but xb is closer to boundries, so xn ends up with being out of boundaries, 
                //so we have to check.
                if (xn <= _lower)
                {
                    xn = _lower;
                }
                if (xn >= _upper)
                {
                    xn = _upper;
                }
                fn = _logfunc(xn);

                //we should never accept fn in here in this 1-d case, since fn is never between fs and fb, since fs and fb in this case is the same.


                //check what to do next
                if (fn > fb)  //so this mean we need to contract, in this case fn<fsecond for other cases, now we need to contract
                { //contract inside or outside???
                    double xc, fc;
                    if (fn < fw) //contract outside. this is better
                    {
                        xc = xCentroid + C_beta * (xn - xCentroid);
                        fc = _logfunc(xc);
                        if (fc < fn)
                        {
                            fn = fc;
                            xn = xc;
                            continue;
                        }
                        //esle will go to do shrink
                    }
                    else //contract inside????
                    {
                        xc = xCentroid + C_beta * (xw - xCentroid);
                        fc = _logfunc(xc);
                        if (fc < fn)
                        {
                            fn = fc;
                            xn = xc;
                            continue;
                        }
                        //else will go to do shrink

                    }
                }
                else //this means fn<fb, so we are doing expansion
                {
                    //compute the expansion point now,
                    double xe, fe;
                    xe = xCentroid + C_gamma * (xn - xCentroid);
                    fe = _logfunc(xe);
                    if (fe <= fn)//will do the expansion
                    {
                        //but need to check whether we are out boundaries,
                        if (xe < _lower)
                        {
                            xe = _lower;
                            fe = _logfunc(xe);
                            if (fe <= fn)//do the expansion
                            {
                                fn = fe;
                                xn = xe;
                                continue;
                            }
                            else //accept the fn xn as the reflection point,
                            {
                                continue;
                            }
                        }
                        if (xe > _upper)
                        {
                            xe = _upper;
                            fe = _logfunc(xe);
                            if (fe <= fn)//do the expansion
                            {
                                fn = fe;
                                xn = xe;
                                continue;
                            }
                            else //accept the fn xn as the reflectio points
                                continue;
                        }
                        fn = fe;
                        xn = xe;
                        continue;
                    }
                    else //fn=fn, fn=xr; accept the reflection point
                    {
                        continue;
                    }
                }

                //if we are here, it means we need to run shrinkage.
                xn = xb + C_sigma * (xw - xb);
                fn = _logfunc(xn);

            }//end of outer while

            return xb;
        }*/

        /*public static double FixBigValueWithFixedIterations(LogDistributionFuctionDelegate _logfunc, double _x0, double _x1, double _lower, double _upper, int _maxIters)
        {
            double f0 = _logfunc(_x0);

            double f1 = _logfunc(_x1);
            NelderMead1D.C_MaximumIterations = _maxIters;
            double xb = FindNonZeroValueKernelLarger(_logfunc, _x0, _x1, _lower, _upper, f0, f0, f1);

            NelderMead1D.C_MaximumIterations=-1;
            return xb;
        }*/

        //declaration of paramter
        static double C_alpha = 1;  //reflection
        static double C_beta = 0.5; //contraction
        static double C_gamma = 2;  //expansion
        static double C_sigma = 0.5; //shrinkage
        static int C_MaximumIterations = 300;//for a search with fixed iterations.
    }//end of class
}//end of namespace
