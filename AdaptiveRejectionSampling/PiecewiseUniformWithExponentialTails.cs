using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Meta.Numerics.Statistics.Distributions;
using System.Numerics;
using NelderMeadMethod;

namespace AdaptiveRejectionSampling
{
    
     
    /// <summary>
    /// this is a simpler way to construct proposal density function based on the paper
    /// 	arXiv:1205.5494 [stat.CO]
    /// 	improved adaptive rejction metrooplis sampling algorithms. 2012. Luca Martina, Jesse Read, David luengo.
    /// this method using uniform way with exponential tail on ends to build the function. it is not assume logConcave and doesn't
    /// lead to "envelope function". so it has to be combined with Adaptive Rejction Metropolis Sampling algorithm (ARMS).
    /// </summary>
    public class PiecewiseUniformWithExponentialTails
    {
        const double ZERO = 1E-20;
        /// <summary>
        /// empty constructor is not allowed for calling (by protected).
        /// </summary>
        protected PiecewiseUniformWithExponentialTails()
        {
            this.InitializeElements();

        }
        /// <summary>
        /// Constructor with support points decided.
        /// </summary>
        /// <param name="_supportPoints">the support points vector, this is on the x axis, so this the abscissae, 
        /// it could be bounded on both end or otherwise. 
        /// in the latter case, it could be minus/plus infinity. In the infinity case, the left end has to be positive gradient
        /// and the right end has to be negative gradient. If the end is bound, then, it is not necessary to have the gradient defined like above.
        /// </param>
        /// <param name="_distributionDesity">the function delegate to the target density distribution</param>
        protected PiecewiseUniformWithExponentialTails(List<double> _supportPoints, LogDistributionFuctionDelegate _logDistributionDesity)
        {
            
            /*CP_SupportPoints = _supportPoints;
            CP_EnvelopeFunctionOfLogTarget = new List<List<double>>();*/

            this.InitializeElements();
            CP_LogDistributionDesity = _logDistributionDesity;
            
            //GenerateEnevelopeFunctionOfLogTarget();
            CreatePiecewiseEnvelopeDistribution(_supportPoints);
        }
        /// <summary>
        /// constructor with number of support points and upper/lower bound specified
        /// </summary>
        /// <param name="_numberOfSupportPoints">number of supporting points/abscissae</param>
        /// <param name="_lowerBound">the lower bound of the support points, could be negative infinity</param>
        /// <param name="_upperBound">the upper bound of the support points, could be positive infinity</param>
        /// <param name="_distributionDesity">the function delegate to the target density distribution</param>
        public PiecewiseUniformWithExponentialTails(int _numberOfSupportPoints, double _lowerBound, double _upperBound, 
                LogDistributionFuctionDelegate _logDistributionDesity, Double _initialValuePrevious, Double _initialValueCurrent)
        {
            //rng1 = new Random();
            /*
            CP_SupportPoints  = new List<double>();
            CP_EnvelopeFunctionOfLogTarget = new List<List<double>>();
            */
            this.C_InitialValue0 = _initialValuePrevious;
            this.C_InitialValue1 = _initialValueCurrent;
            InitializeElements();

            CP_LogDistributionDesity = _logDistributionDesity;

            /*
            //now we need to determine the supportpoints vector.
            GenerateSupportPoints(_numberOfSupportPoints, _lowerBound, _upperBound);
            //now we need to compute the lines and piecewise probability
            GenerateEnevelopeFunctionOfLogTarget();
             */
            CreatePiecewiseEnvelopeDistribution(_numberOfSupportPoints, _lowerBound, _upperBound, _initialValuePrevious,_initialValueCurrent);
        }
        /// <summary>
        /// Generate the supportpoints array based on the number and upper and lower bound
        /// </summary>
        /// <param name="_numberOfSupportPoints">number of abscissae including lower and upper bound, 
        /// so number must larger than 6(including lower and upper bound, or larger than 4 not including lower and upper bound
        /// it is possible to have 5, but it is a little complicate, so we pass this for now.</param>
        /// <param name="_lowerBound">lower bound, could be infinity</param>
        /// <param name="_upperBound">upper bound, could be infinity</param>
        private void GenerateSupportPoints(int _numberOfSupportPoints, double _lowerBound, double _upperBound, double _initial0, double _initial1)
        {
            double originalLowerBound = _lowerBound;
            double originalUpperBound = _upperBound;

            if(_numberOfSupportPoints <=5)
            {
                throw new System.Exception("too few points specified, must be larger then 2");
                //System.Environment.Exit(-1);
            }
            //the hardest point is to make sure, under the infinity bound cases, the two end exponential tails are positive and negative gradient-ed. 
            //for non-infinite bound cases, we don't care. Upated 2013/03/16, actually, we will take care of this case of not positive slope (lower end)
            //or negative slope (upper end) in the createEnvelopefunction method.Here, we simply go ahead, and only search for bigger value on two end than 
            //the boundary points when there is infinity values.
            /*if (Double.IsInfinity(_lowerBound)||Double.IsInfinity(_upperBound))
            {*/
                if(_lowerBound>=_upperBound)
                {
                    throw new System.Exception("the lower and upper bound is not set correctly, lower:"+_lowerBound+" and upper:"+_upperBound );
                }
                //get the smallest possible x1 that makes the function not zero, since zero will make log invalid.
                

            //add new code to take care of the cases where the nonzero region for Target function is too small, which will leave the support array is too small to partition
                while (true)
                {
                    this.CP_SupportPoints = new List<double>();   //we add this to make sure the suppor array now is brand new. this is because sometimes we need to regenerate the support array when we did not found the global optimum the first time which leave us with an overflow of the double values (infinity)


                    double starting;

                    starting = SearchNonZeroPointNM(CP_LogDistributionDesity, _lowerBound, _upperBound, _initial0, _initial1);
                    List<double> xArr = SearchLowerBound(CP_LogDistributionDesity, starting, _lowerBound);
                    double x1 = xArr[0];
                    if (x1 > _upperBound)
                    {
                        throw new System.Exception("not correct function, the smallest possible x value that is not zero is larger than the upperbound");
                    }
                    _lowerBound = xArr[1];
                    if (!double.IsNaN(_lowerBound))//_lowerBound is not an infinity, since we have search for the lower bound in search lower bound
                    {
                        CP_SupportPoints.Add(_lowerBound);
                        _numberOfSupportPoints--;
                    }
                    /*8888888888
                    double lowerFunValue=Math.Exp(CP_LogDistributionDesity(_lowerBound));
                
                    //in this case, the boundary is infinity or the boundary is not infinity, but the funcValue is not usable (0 or NaN), 
                    //so x1 and lower is not the same. We will add a different x1 other than lower.
                    if (Double.IsInfinity(_lowerBound) || ((!Double.IsInfinity(_lowerBound))&&(lowerFunValue  == 0||Double.IsNaN(lowerFunValue )))) 
                        //this is the case where x1 is the _lowerBound, so won't have to add this
                    {
                        */
                    SupportPoints.Add(x1);
                    _numberOfSupportPoints--;
                    _lowerBound = x1;
                    //flagSameLowerBound = false;
                    /*88888}*/
                    double step = 1.01; double x2;
                    if (Double.IsInfinity(_lowerBound))//this is not possible??? since we have use searchLowerBound to make _lowerBound not an infinity
                        //but fine, we will leave it for here for now.
                    {
                        //now we need this to exponential tail to be positive gradient line L(1,2){x;s}
                        //double x1 = -1E1;

                        if (x1 > 0)
                            x2 = x1 * step;
                        else
                        {
                            step = 0.99;
                            x2 = x1 * 0.01;
                        }
                        if (x1 == 0)
                        {
                            x2 = (x1 + _upperBound) / 1000;

                            //need to take care of case, where there is infinity, 
                            //here only upper bound could infinity, since lower can not be infinity for now

                            /*if (Double.IsInfinity(_lowerBound))
                                    x2 = (x1 + _upperBound) * 0.5;
                             */
                            if (Double.IsInfinity(_upperBound))
                                x2 = 0.001;
                        }
                        x2 = NelderMead1D.FindBigValue(this.CP_LogDistributionDesity, _lowerBound, x2, _lowerBound, _upperBound, this.C_FunctionNormConstant);
                        /*        
                        //then check for x2 to get positive slope line
                                while (CP_LogDistributionDesity(x2) < CP_LogDistributionDesity(x1))
                                {
                                    x2 = x1 *step;
                                    if (x2 > _upperBound)
                                    {
                                        throw new System.Exception("there is something wrong with the distribution density funciton, we can not find positive slope line on infinite bound set");
                                    }
                                }
                         */
                        //CP_SupportPoints.Add(_lowerBound);
                        //CP_SupportPoints.Add(x1);
                        CP_SupportPoints.Add(x2);
                        _lowerBound = x2;
                        _numberOfSupportPoints--;


                    }
                    else  //else means lower is not infinite, so then we don't care about what to put in as the x2, since the array is bound, we will take care this case in the createEnvelope function next. 
                    {
                        //CP_SupportPoints.Add(_lowerBound);
                        //CP_SupportPoints.Add(x1);
                        _lowerBound = x1;
                        //_numberOfSupportPoints -= 2;
                        //we also need to make check to make sure the slope at this end is of correct direction, 
                        //otherwise, we need to "replace" the lowerBound to be of the closest zero value, otherwise, the number could possible be too big for calculating the non-normalized probality
                        //to do this, we make a very small step and then seen wether we will have a bigger fuc value than this one 
                    }


                    //***********for upperBound

                    starting = 1E15;

                    starting = SearchNonZeroPointNM(CP_LogDistributionDesity, _lowerBound, _upperBound, _initial0, _initial1);
                    xArr = SearchUpperBound(CP_LogDistributionDesity, starting, _upperBound);
                    x1 = xArr[0];
                    if (x1 < _lowerBound)
                    {
                        throw new System.Exception("not correct function, the smallest possible x value that is not zero is larger than the upperbound");
                    }

                    _upperBound = xArr[1];
                    if (!double.IsNaN(_upperBound))
                    {

                        CP_SupportPoints.Add(_upperBound);
                        _numberOfSupportPoints--;
                    }
                    /*8888888
                        double upperFunValue = Math.Exp(CP_LogDistributionDesity(_upperBound));

                       if (Double.IsInfinity(_upperBound) || ((!Double.IsInfinity(_upperBound))&&(upperFunValue  == 0||Double.IsNaN(upperFunValue )))) //this is the case where x1 is the _lowerBound
                        {*/
                    SupportPoints.Add(x1);
                    _numberOfSupportPoints--;
                    _upperBound = x1;
                    //flagSameUpperBound = false;
                    /*8888}*/
                    step = 0.99;
                    if (Double.IsInfinity(_upperBound))
                    { //in this case, we have to add the point to make the slope of exponential line to be negative

                        if (x1 > 0)
                            x2 = x1 * step;
                        else
                        {
                            step = 1.01;
                            x2 = x1 * step;
                        }
                        if (x1 == 0)
                        {
                            x2 = 0.001 * _lowerBound;
                            //at this point _lowerBound can not be infinity, _lwerBound is the x2 from last round

                        }
                        //find a bigger value on the left side of x1 (the boundary value) to make the slope negative on the boundary
                        x2 = NelderMead1D.FindBigValue(this.CP_LogDistributionDesity, x1, x2, _lowerBound, _upperBound, this.C_FunctionNormConstant);
                        /*
                                //then check for x2 to get positive slope line
                                //bool flag=false;
                                while (CP_LogDistributionDesity(x2) < CP_LogDistributionDesity(x1))
                                {
                                    x2 = x1 * step;
                                    if (x2 < _lowerBound)
                                    {
                                        throw new System.Exception("there is something wrong with the distribution density funciton, we can not find positive slope line on infinite bound set");
                                    }
                                }
                         */
                        //CP_SupportPoints.Add(_upperBound);
                        //CP_SupportPoints.Add(x1);
                        CP_SupportPoints.Add(x2);
                        _upperBound = x2;
                        _numberOfSupportPoints--;
                    }
                    else //in here, the upperBound is not infinity, so, we don't care about what to put in as x2, since it is bounded, we will take care of not create slope in the createEnvelope function next
                    {
                        //CP_SupportPoints.Add(_upperBound);
                        //CP_SupportPoints.Add(x1);
                        //_numberOfSupportPoints -= 2; 
                        _upperBound = x1;
                    }

                    //here we need to check where the support array is too narrow, which means the array is a point.
                    if (
                        (_upperBound == 0 && _lowerBound == 0) ||
                        (_upperBound != 0 && (((_upperBound - _lowerBound) / _upperBound < 1E-7) && ((_upperBound - _lowerBound) / _upperBound > -1E-7))) ||
                        (_lowerBound != 0 && (((_upperBound - _lowerBound) / _lowerBound < 1E-7) && ((_upperBound - _lowerBound) / _lowerBound > -1E-7)))

                    )
                    {
                        this.C_FunctionNormConstant += 10; //relax the nozero region
                        _upperBound = originalUpperBound;
                        _lowerBound = originalLowerBound;
                    }
                    else
                    {
                        break;
                    }
                }//end of while loop to make sure the support array is not too narrow 

                //now we need to add the middle ones.
                for (int i = 1; i <= _numberOfSupportPoints; i++)
                {
                    CP_SupportPoints.Add(_lowerBound + i * (_upperBound - _lowerBound) / (_numberOfSupportPoints + 1));
                }
                CP_SupportPoints.Sort();
            //}
            /*else  //non-infinite bounds
            {
                CP_SupportPoints.Add(_lowerBound);
                
                for (int i = 1; i < _numberOfSupportPoints-1 ; i++)
                {
                    CP_SupportPoints.Add(_lowerBound + i * (_upperBound - _lowerBound) / (_numberOfSupportPoints-1));
                }
                CP_SupportPoints.Add(_upperBound);
            }*/
        }//end of method

        /// <summary>Deprecated!!!! since Mar/4/2013
        /// this is the function to search nonZero point (f(x)!=0, so, we can use this as starting point to the area where we can use this function, 
        /// if in the area where this function value is zero, we can not use log(function) 
        /// </summary>
        /// <param name="function"></param>
        /// <param name="starting"></param>
        /// <param name="step"></param>
        /// <returns></returns>
        /// 
        /*
        public double SearchNonZeroPoint_Depricated(LogDistributionFuctionDelegate _LogFunction,double _lower, double _upper )
        {
            //Console.WriteLine("searching nonzero points");
            //we normall starting with the middle of defined region, or zero if (-inf, +inf), or the one bound if half bound
            double starting = 0;
            double step = 0.001;
            if (!Double.IsInfinity(_lower))
            {
                starting=_lower;
            }
            if (!Double.IsInfinity(_upper))
            {
                starting = _upper;
            }
            if (!Double.IsInfinity(_lower) && !Double.IsInfinity(_upper))
            {
                starting = 0.5 * (_upper + _lower);
                step = (_upper - _lower) / 1E-4;
            }
            //start searching.
            bool found=true;
            int count = 0;
            double originalStarting=starting;

            //we first searching up
            while (Math.Exp(_LogFunction(starting))==0)
            {
                count++;
                //searching up
                starting += step*count;
                if (starting > _upper)
                {
                    found = false;
                    break;
                }
                if (Double.IsInfinity(_upper))
                {//we will stop until some large number
                    if (Double.IsInfinity(_lower))
                    {
                        if (starting > 1E10)
                        {
                            found = false;
                            break;
                        }
                    }
                    else  //lower is not infinity
                    {
                        if ((_lower>0&&starting > _lower + 1E10)||(_lower<0&&starting>1E10))
                        {
                            found = false;
                            break;
                        }
                    }
                }
            }//end of while

            //now we search down

            if (found == false)//we haven't found it through previous searching steps
            {
                starting = originalStarting;
                count = 0;
                found = true;
                while (Math.Exp(_LogFunction(starting))==0)
                {
                    count++;
                    //searching down
                    starting -= step * count;
                    if (starting < _lower)
                    {
                        found = false;
                        break;
                    }
                    if (Double.IsInfinity(_lower))
                    {//we will stop until some small number
                        if (Double.IsInfinity(_upper))
                        {
                            if (starting < -1E10)
                            {
                                found = false;
                                break;
                            }
                        }
                        else  //lower is not infinity
                        {
                            if ((_upper<0)&&(starting < _upper - 1E10)||(_upper>0)&&(starting<-1E10))
                            {
                                found = false;
                                break;
                            }
                        }
                    }

                }//end of the second while
            }//end of second outer if loop
            //now we need to check for whether we found it or not.
            if (found == false)
            {
                //no!!!
                throw new System.Exception("can not found any value on the desity to be nozero");
            }
            Console.WriteLine("done for searching for nonzero.......");
            return starting;
        }//end of function
        */
        /// <summary>
        /// this is a new version of serch non zero point, based on reference
        /// Rajeeva L. Karandikar. On Adaptive Rejction Sampling. 2005. http://www.isid.ac.in/~statmath/eprints
        /// 
        /// NOte: 3/14/2013, this still not working well, so make changes on 3/14/2013 to search around the current value of the function
        /// assume the non-zero value is not far from the current value.
        /// </summary>
        /// <param name="_function">the target function</param>
        /// <returns></returns>
        /*public double SearchNonZeroPointNew_Depricated(LogDistributionFuctionDelegate _LogFunction, double _lower, double _upper, double _initial)
        {
            if (Math.Exp(_LogFunction(_initial)) > 0)
                return _initial;
            //int n = 0;
            for (int n = 1; n<12 ; n++)
            {
                Console.WriteLine("\tsearching.N is " + n);
                double low=intPower(2,4*n);
                //double up=intPower(2, 5*n+1);
                double  step=intPower(2,n);
                for (int i = 1; i*step <= low; i++)
                {
                    double T = _initial + i / step;
                    if(Math.Exp(_LogFunction(T))>0)
                    {
                        Console.WriteLine("\tdone searching!");
                        return T;
                    }
                    T = _initial - i / step;
                    if (Math.Exp(_LogFunction(T)) > 0)
                    {
                        Console.WriteLine("\tdone searching!");
                        return T;
                    }
                }
            }
            return Double.NaN;
        }*/

        /// <summary>
        /// search non zero point using Nelder Mead method
        /// </summary>
        /// <param name="_LogFunction"></param>
        /// <param name="_lower"></param>
        /// <param name="_upper"></param>
        /// <param name="_initial0"></param>
        /// <param name="_initial1"></param>
        /// <returns></returns>
        public double SearchNonZeroPointNM(LogDistributionFuctionDelegate   _LogFunction, double _lower, double _upper, 
            double _initial0, double _initial1)
        {
            return NelderMead1D.FindNonZeroValue(_LogFunction, _initial0, _initial1, _lower, _upper, NelderMead1D.LOG_LIMIT, ref this.C_FunctionNormConstant  );
        }
        /*private double intPower(int _base, int _power)
        {
            double ret = 1.0;
            for (int i = 0; i < _power; i++)
            {
                ret *= _base;
            }
            return ret;
        }*/

        /// <summary>
        /// please set the starting point to be of nonzero value, suppose we started in the middle of the function (nonzero f value), then using 
        /// binary search method. search for the extreme values (or closest to the extreme values) that not equal to zero
        /// </summary>
        /// <param name="function"></param>
        /// <param name="startingPoint">the starting points should be of nonzero f value, </param>
        /// 
        /// <returns>the smallest x that the function will not be zero at position 0, but the one biggest x that will be zero at postion 1
        /// will return nan at position 1 when the lower bound itself is not zero.</returns>
        private List<double> SearchLowerBound(LogDistributionFuctionDelegate _LogFunction, double startingPoint, double _lower)
        {
            List<double> retList = new List<double>(2);
            double f=Math.Exp(_LogFunction(_lower, this.C_FunctionNormConstant ));
            if (f != 0 && !Double.IsNaN(f) && !double.IsInfinity(f))
            {

                retList.Add( _lower);
                retList.Add(double.NaN);
                return retList;
            }
            double step;
            double gap=1E10;
            //double last = f;
            if (Double.IsInfinity(_lower))
            {
                _lower = -1E20;
            }
            while (gap >1E-6)
            {
                step = (startingPoint - _lower) / 2; double temp = Math.Exp(_LogFunction(step + _lower, this.C_FunctionNormConstant ));
                if (temp <= Math.Exp(NelderMead1D.LOG_LIMIT) || Double.IsNaN(temp))
                {
                    _lower=step+_lower;
                }
                else
                {
                    startingPoint=step+_lower;
                }
                gap = (startingPoint - _lower)/startingPoint;
            }
            retList.Add(startingPoint);
            retList.Add(_lower);
            return retList ;
        }
        
        /// <summary>
        /// please set this starting point to be nonzero value, suppose we are starting with nonzero f value starting point, then use binary search to get the
        /// largest or close to larget nonzero f value abscissae
        /// </summary>
        /// <param name="function"></param>
        /// <param name="startingPoint">non zero f value x abscissae</param>
        /// 
        /// <returns>at position 0, the biggest value that is not zero;position 1, smallest value that is zero. or NaN when the both are upper, or the upper is not zero itself.</returns>
        private List<double> SearchUpperBound(LogDistributionFuctionDelegate _LogFunction, double startingPoint, double _upper)
        {
            List<double> retList = new List<double>(2);
            double f =Math.Exp( _LogFunction(_upper, this.C_FunctionNormConstant));
            if (f != 0 && !double.IsNaN(f) && !double.IsInfinity(f))
            {
                retList.Add(_upper);
                retList.Add(double.NaN);
                return retList;
            }
            double step;
            double gap = 1E10;
            if (Double.IsInfinity(_upper))
            {
                _upper = 1E20;
            }
            //double last = f;
            while (gap > 1E-6||gap<-1E-6)
            {
                step = (  _upper -startingPoint) / 2;
                double temp =Math.Exp( _LogFunction(step + startingPoint,C_FunctionNormConstant ));
                if ( temp<=Math.Exp(NelderMead1D.LOG_LIMIT)||Double.IsNaN(temp))
                {
                    _upper = step+startingPoint;
                }
                else
                {
                    startingPoint = step+startingPoint;
                }
                gap =( _upper- startingPoint )/startingPoint ;
            }
            retList.Add(startingPoint);
            retList.Add(_upper);
            return retList ;
        }
        /// <summary>
        /// this is the function called after a new point is inserted in the previously existing support array, 
        /// we have now known the _index where we put the new point and also we want to call to update only the 
        /// affected points, not necessarily call all other points.
        /// </summary>
        /// <param name="_index">the index where the new point inserted in the support array, it can not be first or last position</param>
        private void UpdateEnevelopeFunctionOfLogTarget(int _index)
        {
            //if the support points array is of N size, this enevelope function is of size N-1
            //so far the suport point array has been updated (the size is enlarged by 1, but the envelop function has not been updated

            //if there is zero content, something wrong, it should have some content.
            if (this.CP_EnvelopeFunctionOfLogTarget==null||CP_EnvelopeFunctionOfLogTarget.Count == 0)
            {
                throw new System.Exception("something wrong here, the envelop function of logTarget has not been ininilized correctly");

                //CP_EnvelopeFunctionOfLogTarget = new List<List<double>>();
            }
            
            //add one more point to the end,************* 
            double logTargetValueOfNewPoint=this.CP_LogDistributionDesity(this.CP_SupportPoints[_index], this.C_FunctionNormConstant);
            List<double> logTargetEnvelopLine=new List<double>(2);
            logTargetEnvelopLine.Add(0);
            logTargetEnvelopLine.Add(logTargetValueOfNewPoint);//so far assume a horizontal line
            this.CP_EnvelopeFunctionOfLogTarget.Add(logTargetEnvelopLine );
            
            
            //assume the support points is [s0, s1, s2, ....., sn-1, sn, sn+1, NewLine] (N'=n+2 + 1)
            //int i = 1;
            double x1, y1, x2, y2;

            int indexToBeShifted = 0; List<double> OriginalValueInPreviousEnvelopFuctionLine = new List<double>(2) ;
            //*********LEFT TAIL most***************
            if (_index == 1)
            {
                //first, we need to add the left tail of exponential line(the line before exponential, the s1 and s2 line for s0<x<s1

                y2 = CP_LogDistributionDesity(CP_SupportPoints[_index + 1], this.C_FunctionNormConstant);
                y1 = logTargetValueOfNewPoint;//CP_LogDistributionDesity(CP_SupportPoints[i], this.C_FunctionNormConstant);
                x2 = CP_SupportPoints[_index + 1];
                x1 = CP_SupportPoints[_index];

                double slope = (y2 - y1) / (x2 - x1);


                double intercept = y2 - x2 / (x1 - x2) * (y1 - y2);


                List<double> line = new List<double>(2);
                line.Add(slope); line.Add(intercept);
                CP_EnvelopeFunctionOfLogTarget[0]=line;

                //now we need to update the envelop point [1]
                double temp3rdPointLogTarget = y2;//this.CP_LogDistributionDesity(this.CP_SupportPoints[2], this.C_FunctionNormConstant);
                indexToBeShifted = 2;
                OriginalValueInPreviousEnvelopFuctionLine = this.CP_EnvelopeFunctionOfLogTarget[1];
                //List<double> temp2ndEnvelopeFunctionLine = this.CP_EnvelopeFunctionOfLogTarget[1];
                if (logTargetValueOfNewPoint > temp3rdPointLogTarget )
                {
                    CP_EnvelopeFunctionOfLogTarget[1] = logTargetEnvelopLine;
                }
                else
                {
                    CP_EnvelopeFunctionOfLogTarget[1] = new List<double>() { 0, temp3rdPointLogTarget };
                }
                /*******88888888888888888888888888888*****************
                //now we need to shift towards the end by copy over the one before it.
                for (int j = 2; j < this.CP_EnvelopeFunctionOfLogTarget.Count; j++)
                {
                    List<double> tempTemp = this.CP_EnvelopeFunctionOfLogTarget[j];
                    this.CP_EnvelopeFunctionOfLogTarget[j] = temp2ndEnvelopeFunctionLine;
                    temp2ndEnvelopeFunctionLine = tempTemp;
                }*/
            }

            //************left tail second
            if (_index == 2)
            {
                //first, we need to add the left tail of exponential line(the line before exponential, the s1 and s2 line for s0<x<s1

                y2 = logTargetValueOfNewPoint;
                y1 = CP_LogDistributionDesity(CP_SupportPoints[_index - 1], this.C_FunctionNormConstant);//CP_LogDistributionDesity(CP_SupportPoints[i], this.C_FunctionNormConstant);
                x2 = CP_SupportPoints[_index ];
                x1 = CP_SupportPoints[_index-1];

                double slope = (y2 - y1) / (x2 - x1);


                double intercept = y2 - x2 / (x1 - x2) * (y1 - y2);


                List<double> line = new List<double>(2);
                line.Add(slope); line.Add(intercept);
                CP_EnvelopeFunctionOfLogTarget[0]=line;

                //now we need to update the envelop point [1]
                double tempPointLogTarget = y1;//this.CP_LogDistributionDesity(this.CP_SupportPoints[_index-1], this.C_FunctionNormConstant);
                if (logTargetValueOfNewPoint > tempPointLogTarget )
                {
                    CP_EnvelopeFunctionOfLogTarget[1] = logTargetEnvelopLine;
                }
                else
                {
                    CP_EnvelopeFunctionOfLogTarget[1] = new List<double>() { 0, tempPointLogTarget };
                }
                //now we need to update the one after it
                indexToBeShifted = 3;
                OriginalValueInPreviousEnvelopFuctionLine = this.CP_EnvelopeFunctionOfLogTarget[2];
                tempPointLogTarget = this.CP_LogDistributionDesity(this.CP_SupportPoints[3], this.C_FunctionNormConstant);
                if (logTargetValueOfNewPoint > tempPointLogTarget )
                {
                    CP_EnvelopeFunctionOfLogTarget[2] = logTargetEnvelopLine;
                }
                else
                {
                    CP_EnvelopeFunctionOfLogTarget[2] = new List<double>() { 0, tempPointLogTarget };
                }
                //List<double> temp2ndEnvelopeFunctionLine = this.CP_EnvelopeFunctionOfLogTarget[1];
                
                
            }

            //*************in the middle
            //for the case, the point is inserted in the middle
            if (_index > 2 && _index < this.CP_SupportPoints.Count - 3)
            {
                //we need to compare the one before it and the one after it
                //first, compare the one before it
                double tempPointLogTarget = this.CP_LogDistributionDesity(this.CP_SupportPoints[_index-1], this.C_FunctionNormConstant);
                List<double> tempEnvelopeFunctionLine = new List<double>(2);
                if (logTargetValueOfNewPoint > tempPointLogTarget )
                {
                    tempEnvelopeFunctionLine = this.CP_EnvelopeFunctionOfLogTarget[_index - 1];
                    this.CP_EnvelopeFunctionOfLogTarget[_index - 1] = logTargetEnvelopLine ;
                }
                else  // update the one before it 
                {
                    tempEnvelopeFunctionLine = logTargetEnvelopLine;//it won't be used anyway
                    this.CP_EnvelopeFunctionOfLogTarget[_index - 1] = new List<double>() { 0, tempPointLogTarget };
                }

                //to compare with the point AFTER it
                tempPointLogTarget = this.CP_LogDistributionDesity(this.CP_SupportPoints[_index + 1], this.C_FunctionNormConstant);
                tempEnvelopeFunctionLine = this.CP_EnvelopeFunctionOfLogTarget[_index];

                indexToBeShifted = _index+1;
                OriginalValueInPreviousEnvelopFuctionLine = this.CP_EnvelopeFunctionOfLogTarget[_index];
                if (logTargetValueOfNewPoint > tempPointLogTarget )
                {
                    //tempEnvelopeFunctionLine = this.CP_EnvelopeFunctionOfLogTarget[_index ];
                    this.CP_EnvelopeFunctionOfLogTarget[_index ] = logTargetEnvelopLine;
                }
                else  //update with the point after it 
                {
                    //tempEnvelopeFunctionLine = logTargetEnvelopLine;//it won't be used anyway
                    this.CP_EnvelopeFunctionOfLogTarget[_index ] = new List<double>() { 0, tempPointLogTarget };
                }

                //now update the resest of it.
            }

            //*************right tail second to the last
            if (_index == this.CP_SupportPoints.Count - 3)
            {
                //compare with the one before the current point in the support array
                double tempPointLogTarget = this.CP_LogDistributionDesity(this.CP_SupportPoints[_index - 1], this.C_FunctionNormConstant);
                List<double> tempEnvelopeFunctionLine = new List<double>(2);
                if (logTargetValueOfNewPoint > tempPointLogTarget)
                {
                    //tempEnvelopeFunctionLine = this.CP_EnvelopeFunctionOfLogTarget[_index - 1];
                    this.CP_EnvelopeFunctionOfLogTarget[_index-1 ] = logTargetEnvelopLine;
                }
                else  // update the one before it 
                {
                    //tempEnvelopeFunctionLine = logTargetEnvelopLine;//it won't be used anyway
                    this.CP_EnvelopeFunctionOfLogTarget[_index-1 ] = new List<double>() { 0, tempPointLogTarget };
                }

                //compare with the one after it
                tempPointLogTarget = this.CP_LogDistributionDesity(this.CP_SupportPoints[_index + 1], this.C_FunctionNormConstant);
                //List<double> tempEnvelopeFunctionLine = new List<double>(2);
                if (logTargetValueOfNewPoint > tempPointLogTarget)
                {
                    //tempEnvelopeFunctionLine = this.CP_EnvelopeFunctionOfLogTarget[_index - 1];
                    this.CP_EnvelopeFunctionOfLogTarget[_index] = logTargetEnvelopLine;
                }
                else  // update the one before it 
                {
                    //tempEnvelopeFunctionLine = logTargetEnvelopLine;//it won't be used anyway
                    this.CP_EnvelopeFunctionOfLogTarget[_index] = new List<double>() { 0, tempPointLogTarget };
                }

                //doing the tail line
                y2 = tempPointLogTarget;//CP_LogDistributionDesity(CP_SupportPoints[_index+1], this.C_FunctionNormConstant);
                y1 = logTargetValueOfNewPoint;//CP_LogDistributionDesity(CP_SupportPoints[i], this.C_FunctionNormConstant);
                x2 = CP_SupportPoints[_index+1];
                x1 = CP_SupportPoints[_index ];

                double slope = (y2 - y1) / (x2 - x1);


                double intercept = y2 - x2 / (x1 - x2) * (y1 - y2);


                List<double> line = new List<double>(2);
                line.Add(slope); line.Add(intercept);
                CP_EnvelopeFunctionOfLogTarget[_index+1] = line;

                indexToBeShifted = _index +2;
                OriginalValueInPreviousEnvelopFuctionLine = null;//at this point, we already reach the end of the array, so we won't update this anyway.

            }

            //*************right tail most
            if (_index == this.CP_SupportPoints.Count -2)
            {
                //compare with the one before the current point in the support array
                double tempPointLogTarget = this.CP_LogDistributionDesity(this.CP_SupportPoints[_index - 1], this.C_FunctionNormConstant);
                List<double> tempEnvelopeFunctionLine = new List<double>(2);
                if (logTargetValueOfNewPoint > tempPointLogTarget)
                {
                    //tempEnvelopeFunctionLine = this.CP_EnvelopeFunctionOfLogTarget[_index - 1];
                    this.CP_EnvelopeFunctionOfLogTarget[_index-1] = logTargetEnvelopLine;
                }
                else  // update the one before it 
                {
                    //tempEnvelopeFunctionLine = logTargetEnvelopLine;//it won't be used anyway
                    this.CP_EnvelopeFunctionOfLogTarget[_index-1] = new List<double>() { 0, tempPointLogTarget };
                }

                

                //doing the tail line
                y2 = tempPointLogTarget;//CP_LogDistributionDesity(CP_SupportPoints[_index - 1], this.C_FunctionNormConstant);
                y1 = logTargetValueOfNewPoint;//CP_LogDistributionDesity(CP_SupportPoints[i], this.C_FunctionNormConstant);
                x2 = CP_SupportPoints[_index -1 ];
                x1 = CP_SupportPoints[_index];

                double slope = (y2 - y1) / (x2 - x1);


                double intercept = y2 - x2 / (x1 - x2) * (y1 - y2);


                List<double> line = new List<double>(2);
                line.Add(slope); line.Add(intercept);
                CP_EnvelopeFunctionOfLogTarget[_index ] = line;

                indexToBeShifted = _index + 1;
                OriginalValueInPreviousEnvelopFuctionLine = null;//at this point, we already reach the end of the array, so we won't update this anyway.


            }
            //now we need to shift the envelope array to the end if there are any.
            //*******88888888888888888888888888888*****************
            //now we need to shift towards the end by copy over the one before it.
            for (int k = indexToBeShifted; k < this.CP_EnvelopeFunctionOfLogTarget.Count; k++)
            {
                List<double> tempTemp = this.CP_EnvelopeFunctionOfLogTarget[k];
                this.CP_EnvelopeFunctionOfLogTarget[k] = OriginalValueInPreviousEnvelopFuctionLine ;
                OriginalValueInPreviousEnvelopFuctionLine  = tempTemp;
            }
            /*for (int j = 1; j < CP_SupportPoints.Count - 2; j++)
            {
                line = new List<double>(2);
                line.Add(0);//always 
                line.Add(
                    Math.Max(CP_LogDistributionDesity(CP_SupportPoints[j], this.C_FunctionNormConstant), CP_LogDistributionDesity(CP_SupportPoints[j + 1], this.C_FunctionNormConstant))
                    );
                CP_EnvelopeFunctionOfLogTarget.Add(line);
            }

            //now add the right tail; the right exponential tail (the line before exponential, the sn-1, sn line for sn<x<sn+1
            i = CP_SupportPoints.Count - 3;
            x1 = CP_SupportPoints[i];
            x2 = CP_SupportPoints[i + 1];
            y1 = CP_LogDistributionDesity(CP_SupportPoints[i], this.C_FunctionNormConstant);
            y2 = CP_LogDistributionDesity(CP_SupportPoints[i + 1], this.C_FunctionNormConstant);

            slope = (y2 - y1) / (x2 - x1);
            intercept = y2 - x2 / (x1 - x2) * (y1 - y2);
            line = new List<double>(2);
            line.Add(slope); line.Add(intercept);
            CP_EnvelopeFunctionOfLogTarget.Add(line);*/
        }

        private void GenerateEnevelopeFunctionOfLogTarget()
        {//we need to go through the supportPoints array to get the envelopefunction of the log target
            //if the support points array is of N size, this enevelope function is of size N-1

            //if this not zero content, then it means we are doing updating/insertion.
            if (CP_EnvelopeFunctionOfLogTarget.Count > 0)
            {
                CP_EnvelopeFunctionOfLogTarget = new List<List<double>>();
            }
            //assume the support points is [s0, s1, s2, ....., sn-1, sn, sn+1] (N=n+2)
            //*********LEFT TAIL***************
            //first, we need to add the left tail of exponential line(the line before exponential, the s1 and s2 line for s0<x<s1
            int i = 1;
            double x1, y1, x2, y2;
            y2=CP_LogDistributionDesity(CP_SupportPoints[i + 1], this.C_FunctionNormConstant );
            y1=CP_LogDistributionDesity(CP_SupportPoints[i ], this.C_FunctionNormConstant);
            x2=CP_SupportPoints[i + 1] ;
            x1= CP_SupportPoints[i];

            double slope = (y2-y1)  / (x2-x1);
            
            
            double intercept = y2 - x2 / (x1 - x2) * (y1 - y2);
            
                        
            List<double> line = new List<double>(2);
            line.Add(slope); line.Add(intercept);
            CP_EnvelopeFunctionOfLogTarget.Add(line);

            //now we need to add the horizontal straight line
            for (int j = 1; j < CP_SupportPoints.Count - 2; j++)
            {
                line = new List<double>(2);
                line.Add(0);//always 
                line.Add( 
                    Math.Max(CP_LogDistributionDesity(CP_SupportPoints[j], this.C_FunctionNormConstant ), CP_LogDistributionDesity(CP_SupportPoints[j+1], this.C_FunctionNormConstant ))
                    );
                CP_EnvelopeFunctionOfLogTarget.Add(line);
            }

            //now add the right tail; the right exponential tail (the line before exponential, the sn-1, sn line for sn<x<sn+1
            i = CP_SupportPoints.Count - 3;
            x1=CP_SupportPoints[i];
            x2=CP_SupportPoints[i+1];
            y1=CP_LogDistributionDesity(CP_SupportPoints[i], this.C_FunctionNormConstant );
            y2=CP_LogDistributionDesity(CP_SupportPoints[i+1], this.C_FunctionNormConstant );

            slope = (y2-y1) / (x2 - x1);
            intercept = y2 - x2 / (x1 - x2) * (y1 - y2);
            line = new List<double>(2);
            line.Add(slope); line.Add(intercept);
            CP_EnvelopeFunctionOfLogTarget.Add(line);
        }

        /// <summary>
        /// this is generate probablity for each piece based on exp( H function), also nomalized by the total probablity of the sum
        /// Updated on Apr 4th 2013, here we add changes to take care the case where we have overflow (infinity) for piecewise PDF value
        /// so far, the reason for this is that we haven't found the correct maximum value of the curve (log target function), so the 
        /// NelderMead could be stuck in the local optimum, that is why. if we are luck we could include the "global" optimum in the 
        /// support array, and then we are out of infinity when dealing with this optimal points. To solve this problem, we can check for the new optimum 
        /// using the biggest value found in this array. this will give us better one or hopefully the best.
        /// </summary>
        /// <param name="_H">envelope function, H</param>
        /// <param name="_S">support points array</param>
        /// <returns></returns>
        public List<double> GeneratePiecewiseProbabilityOfEnvelopeFunction(List<List<double>> _H, List<double> _S)
        {

            List<double> p;
            double p_total = 0;
            int indexBiggestLogTargetFunctionValue = 1;//we don't care about the first [0] and last [Count-1] values, since they are only for bound purpose, not having meaningful values for our purpose
            double biggestLogTargetFunctionValue = CP_EnvelopeFunctionOfLogTarget[indexBiggestLogTargetFunctionValue][1];

            while (true)
            {
                p = new List<double>(_H.Count);
                //first, for the left tail, _h : ax+b, for exponential, this is exp(ax+b), where a >0, and x<-(-inf, s1)
                //integration of this : f(exp(ax+b)dx) =>1/a*exp(ax+b)|[s1, -inf], f is the integral  ==>1/a*exp(a*s1-b)
                if (_H[0][0] != 0)
                {
                    p.Add(1 / _H[0][0] * (Math.Exp(_H[0][0] * _S[1] + _H[0][1]) - Math.Exp(_H[0][0] * _S[0] + _H[0][1])));
                }
                else
                {
                    if (Double.IsInfinity(_S[0]))
                        throw new System.Exception("the envelope function has not been constructed correctly, the tail is a flat with a infinite lower bound");
                    p.Add(Math.Exp(_H[0][1]) * (_S[1] - _S[0]));
                }
                p_total = p[0];

                //doing the middle ones, uniform in exp format
                for (int i = 1; i < _H.Count - 1; i++)
                {
                    p.Add(Math.Exp(_H[i][1]) * (_S[i + 1] - _S[i]));
                    p_total += p[i];//the one just added, also could be p[p.count-1]
                    if (_H[i][1] > biggestLogTargetFunctionValue)
                    {
                        indexBiggestLogTargetFunctionValue = i;
                        biggestLogTargetFunctionValue = _H[i][1];
                    }
                }

                //right tail, _h: ax+b, for expontial, this exp(ax+b), where a<0, and x<-(sn, -inf), integration of 
                //this:  f(exp(ax+b)dx) =>1/a*exp(ax+b)|[+inf, Sn], f is the integral  ==>-1/a*exp(a*s1-b)
                if (_H[_H.Count - 1][0] != 0)
                {
                    p.Add(1 / _H[_H.Count - 1][0] * (Math.Exp(_H[_H.Count - 1][0] * _S[_S.Count - 1] + _H[_H.Count - 1][1]) -
                    Math.Exp(_H[_H.Count - 1][0] * _S[_S.Count - 2] + _H[_H.Count - 1][1])));
                }
                else
                {
                    if (Double.IsInfinity(_S[_S.Count - 1]))
                        throw new System.Exception("the envelope function has not been constructed correctly, the tail is a flat with a infinite upper bound");
                    p.Add(Math.Exp(_H[_H.Count - 1][1]) * (_S[_S.Count - 1] - _S[_S.Count - 2]));

                }

                p_total += p[p.Count - 1];
                if (Double.IsInfinity(p_total))
                {
                    this.C_FunctionNormConstant = 0;
                    //now we need to call to generate the support point array and envelope again
                    this.GenerateSupportPoints(this.CP_SupportPoints.Count, this.CP_SupportPoints[0], this.CP_SupportPoints[this.CP_SupportPoints.Count - 1], this.CP_SupportPoints[indexBiggestLogTargetFunctionValue - 1], this.CP_SupportPoints[indexBiggestLogTargetFunctionValue]);
                    this.GenerateEnevelopeFunctionOfLogTarget();
                    _S = this.CP_SupportPoints;
                    _H = this.CP_EnvelopeFunctionOfLogTarget;
                }
                else
                {//we are good so far, so go ahead.
                    break;
                }
            }
            //now we piece normalized probablity, after this it should add up to 1.
            for (int j = 0; j < p.Count; j++)
            {
                p[j] /= p_total;
                
            }
            CP_IntegralOfTotalEnvelopeFunction = p_total;
            return p;
        }
        /// <summary>
        /// this is the evaluated value of exponential envelope function.used to compare with function(x) on ARS
        /// </summary>
        /// <param name="_x"></param>
        /// <returns></returns>
        public double ExpHnOfX(double _x)
        {
            //first we need to find the region where this_x belongs to, then retunr the expHnOfX
            int index=-1;
            for (int i=0;i<CP_SupportPoints.Count-1;i++)
            {
                if (_x >= this.CP_SupportPoints[i] && _x < this.CP_SupportPoints[i + 1])
                {
                    index = i ;
                    break;
                }
            }
            if (index ==-1)//this only means we can not find the point in the support array, so this means the expHnOfX is zero. nothing to worry so far.
            {
                //Be careful here!!!!
                //Console.WriteLine ("WARNING.........x doesn't belong the defined envelope function");
                return 0;
            }
            return Math.Exp(this.CP_EnvelopeFunctionOfLogTarget[index][0] * _x + this.CP_EnvelopeFunctionOfLogTarget[index][1]);
        }

        public List<double> GeneratePiecewiseCDFOfEnvelopeFunction(List<double> _pdf)
        {
            List<double> cdf = new List<double>(_pdf.Count);
            cdf.Add(_pdf[0]);
            double total_p = _pdf[0];
            for (int j = 1; j < _pdf.Count; j++)
            {
                total_p+=_pdf[j];
                cdf.Add(total_p);
            }
            return cdf;
        }


        //now get the sampling function
        public double GetRandomValue(Random _rng)
        {
            UniformDistribution ufd = new UniformDistribution();
            double rn = ufd.GetRandomValue(_rng);

            //Console.WriteLine("random number is " + rn);
            //using the inversion method to generate a sample from the proposal/envelope function
            //first go through the piece wise cdf to get the "piece" first
            int piece = 0;
            for (int i = 0; i < CP_PiecewiseCDF.Count; i++)
            {
                if (rn < CP_PiecewiseCDF[i])
                {
                    piece = i;
                    break;
                }
            }
            if (piece > 0)
                rn -= CP_PiecewiseCDF[piece-1];

            if (piece == 0 || piece == CP_PiecewiseCDF.Count - 1)
            {//this is tail, so use exponential
                //1/a*(exp(a*s_i+1 + b) - exp(a*s_i +b)) / integral
                double temp;
                if (CP_EnvelopeFunctionOfLogTarget[piece][0] != 0)
                {
                    temp = rn * CP_IntegralOfTotalEnvelopeFunction * CP_EnvelopeFunctionOfLogTarget[piece][0] + Math.Exp(CP_EnvelopeFunctionOfLogTarget[piece][0] * CP_SupportPoints[piece] + CP_EnvelopeFunctionOfLogTarget[piece][1]);
                    return (Math.Log(temp) - CP_EnvelopeFunctionOfLogTarget[piece][1]) / CP_EnvelopeFunctionOfLogTarget[piece][0];
                }
                else //for special case where the tails are also the uniform line too
                {
                    return rn * CP_IntegralOfTotalEnvelopeFunction / Math.Exp(CP_EnvelopeFunctionOfLogTarget[piece][1]) + CP_SupportPoints[piece];
                }
                
            }
            else //this is the uniform, so use uniform
            {
                //exp(H_i)(S_i+1 -S_i)/integral
                return rn * CP_IntegralOfTotalEnvelopeFunction / Math.Exp(CP_EnvelopeFunctionOfLogTarget[piece][1]) + CP_SupportPoints[piece]; 
            }

        }


        //now need to update the support points array and envelope array and pdf cdf
        //was before Apr 1st 2013.
        public void AddOnePointToSupportPointsArray_old(double _x)
        {
            //first find the index
            int i = 0;
            for (i = 0; i < CP_SupportPoints.Count - 1; i++)
            {
                //check whether _x is between this current region
                if (_x >= CP_SupportPoints[i] && _x < CP_SupportPoints[i + 1])
                {
                    break;
                }
            }
            if (i == CP_SupportPoints.Count - 1) //this is not right
                throw new System.Exception("the newly added point is not in the range of the support points");

            CP_SupportPoints.Insert(i+1, _x);

            UpdatePiecewiseEnvelopeDistribution();
        }
        //now need to update the support points array and envelope array and pdf cdf
        //updated on April 1st 2013 to make it fast and also calling the downstream updateEnvelopeFunction 
        //by not calling what is not necessary
        public void AddOnePointToSupportPointsArray(double _x)
        {
            //first to see whether this current one is between the first and last member of support array.
            //if not, then something is wrong. it also assumes that the array is sorted correctly
            if (_x >= this.CP_SupportPoints[this.CP_SupportPoints.Count - 1] || _x <= this.CP_SupportPoints[0])
            {
                throw new System.Exception("the newly added point is not in the range of the support points");
            }

            //we are so far
            //need to add this point to the end of the support array first
            CP_SupportPoints.Add(_x);
            
            //now we using the swap method to put it in the current position
            int i = this.CP_SupportPoints.Count-2;//so far the last one is _x itself
            int swapTempIndex=this.CP_SupportPoints.Count-1;
            
            for (; i >=0; i--)
            {
                //check whether _x is between this current region
                if (this.CP_SupportPoints[swapTempIndex] < CP_SupportPoints[i])
                {
                    this.CP_SupportPoints[swapTempIndex] = CP_SupportPoints[i];
                    this.CP_SupportPoints[i] = _x;
                    swapTempIndex = i;

                }
                else //we are done now, because in a sorted array, the current value is not smaller anymore.
                {
                    i++;
                    break;
                }
            }
            if (i == CP_SupportPoints.Count - 1||i==0) //this is not right
                throw new System.Exception("the newly added point is not in the range of the support points");

            //CP_SupportPoints.Insert(i + 1, _x);

            //new added function to update the enevelope function log target
            UpdateEnevelopeFunctionOfLogTarget(i);

            //now call to update the piecewise pdf and cdf
            UpdatePiecewiseEnvelopeDistribution();
        }


        /// <summary>
        /// to initialize all the arrays and other elements too.
        /// </summary>
        private void InitializeElements()
        {
            CP_SupportPoints = new List<double>();
            CP_PiecewiseCDF = new List<double>();
            CP_PiecewisePDF = new List<double>();

            CP_EnvelopeFunctionOfLogTarget = new List<List<double>>();

            //rng1 = new Random();
        }

        /// <summary>
        /// this is the function to call individual functions to make the distribution available to be used for sampling for the FIRST time, 
        /// NOT called to update the distribution after the insertiong of
        /// new point in the support point array. assume the support point array is ready.
        /// </summary>
        /// <param name="_initial">is the intialize value when generate this distribution, this one is used for searching nonzero values</param>
        private void CreatePiecewiseEnvelopeDistribution(int _numberOfSupportPoints, double _lower, double _upper, double _initialPrevious, double _initialCurrent)
        {
            //first we need to generate the support points array
            GenerateSupportPoints(_numberOfSupportPoints, _lower, _upper, _initialPrevious, _initialCurrent  );
            GenerateEnevelopeFunctionOfLogTarget();
            CP_PiecewisePDF= GeneratePiecewiseProbabilityOfEnvelopeFunction(CP_EnvelopeFunctionOfLogTarget, CP_SupportPoints);
            CP_PiecewiseCDF = GeneratePiecewiseCDFOfEnvelopeFunction(CP_PiecewisePDF);

            //it is ready to run sampling
        }

        /// <summary>
        /// overloaded with predefined support point array.
        /// </summary>
        /// <param name="_suppointPoints"></param>
        private void CreatePiecewiseEnvelopeDistribution(List<double> _suppointPoints)
        {
            //we have had the array, 
            //so go ahead to call to create envelope function etc.
            GenerateEnevelopeFunctionOfLogTarget();
            CP_PiecewisePDF = GeneratePiecewiseProbabilityOfEnvelopeFunction(CP_EnvelopeFunctionOfLogTarget, CP_SupportPoints);
            CP_PiecewiseCDF = GeneratePiecewiseCDFOfEnvelopeFunction(CP_PiecewisePDF);
        }


        /// <summary>
        /// called after the insertion of the new point, not for the first time after creating the objects
        /// </summary>
        private void UpdatePiecewiseEnvelopeDistribution()
        {
            //assume the support point array has been called/updated. now we need to recalculate the distribution information and make it ready to run sampling
            

            //updated on Apr 1st 2013,*************** 
            //assume now updating EnvelopeFuncitonOfLogTarget is called in add one point to support array.

            //GenerateEnevelopeFunctionOfLogTarget();
            //**************************
            CP_PiecewisePDF = GeneratePiecewiseProbabilityOfEnvelopeFunction(CP_EnvelopeFunctionOfLogTarget, CP_SupportPoints);
            CP_PiecewiseCDF = GeneratePiecewiseCDFOfEnvelopeFunction(CP_PiecewisePDF);
        }
        //********declaration of member variables
        private LogDistributionFuctionDelegate CP_LogDistributionDesity;
        private List<double> CP_SupportPoints;
        private List<double> CP_PiecewiseCDF;
        private List<double> CP_PiecewisePDF;
        private List<List<double>> CP_EnvelopeFunctionOfLogTarget;
        private double CP_IntegralOfTotalEnvelopeFunction; //this is the normalization factor for the probablity of each piece
        private double C_InitialValue0;
        private double C_InitialValue1;
        private double C_FunctionNormConstant = 0;
       

        //****************PROPERTIES******************
        //properties
        public List<double> SupportPoints
        {
            get { return CP_SupportPoints ; }
            //set { seconds = value * 3600; }
        }
        public double FunctionNormConstant
        {
            get { return C_FunctionNormConstant ; }
            //set { seconds = value * 3600; }
        }
        public List<List<double>> EnvelopeFunctionOfLogTarget
        {
            get { return CP_EnvelopeFunctionOfLogTarget; }
        }

        public double NormalizationFactorOfProb
        {
            get { return CP_IntegralOfTotalEnvelopeFunction; }
        }

        public List<double> PiecewisePDF
        {
            get { return CP_PiecewisePDF; }
        }

        public List<double> PiecewiseCDF
        {
            get { return CP_PiecewiseCDF ; }
        }
    }//end of class
}//end of namespace
