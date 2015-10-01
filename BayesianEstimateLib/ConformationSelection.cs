using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BayesianEstimateLib
{
    /// <summary>
    /// this is the basic conformation selection model, without kM or offset. also in this model, the R0_B, R0_B_star, R0_AB and R0_AB_star are all determined by the end of attaching phase.
    /// it states that there are two different conformations for the free antibody and these two are in a steady states
    /// A+B->AB    : ka1, kd1
    /// A+B*->AB*  : ka2, kd2
    /// B->B*      : ka3, kd3
    /// we have ka1,ka2, kd1,kd2, ka3,kd3
    /// different equations
    ///  dAB/dt=ka1xAxB-kd1xAB 
    ///  dAB*/dt=ka1xAxB*-kd2*AB
    ///  dB/dt=-1x dAB/dt -ka3B+kd3B*
    ///  dB*/dt=-1xdAB*/dt - ka3B-kd3B*
    ///  Bmax=B+B*+AB+AB*
    ///  
    /// </summary>
    public class ConformationSelection:SimulationSPR 
    {
          

        protected double _ka2;
        protected double _kd2;

        protected double _ka3;
        protected double _kd3;

        //protected double R0_AB;  <<==this is using the base class SSPR_r0
        protected double _R0_AB_Star;

        protected double _R0_B;
        protected double _R0_B_Star;

        protected double[] _R_AB_Att;
        protected double[] _R_AB_Star_Att;

        protected double[] _R_B_Att;
        protected double[] _R_B_Star_Att;

        protected double[] _R_AB_Det;
        protected double[] _R_AB_Star_Det;

        protected double[] _R_B_Det;
        protected double[] _R_B_Star_Det;

        /*protected double[] _Association_RUs;
        protected double[] _Association_Times;
        protected double[] _Dissociation_RUs;
        protected double[] _Dissociation_Times;

        protected double _dt; */
        
        /*public SPR_Model_TwoState(double[] _ts_association, double[] _ts_dissociation)
        {
            this._times=_ts;
            _dt = -1;
            //do the 
        }*/
        public ConformationSelection()
            : base()
        {
            //empty constructor with unspecified parameters
            //set default, uninitialized values;
            this._ka = -1;
            this._kd = -1;
            this._ka2 = -1;
            this._kd2 = -1;
            this._ka3=-1;
            this._kd3=-1;
            this._conc = -1;
            this._Rmax = -1;
            this.SSPR_r0 = -1;
            this._R0_AB_Star = -1;
            this._R0_B=-1;
            this._R0_B_Star=-1;
            this._duration_attach = -1;
            this._duration_detach = -1;
            //this._deltaT =;
            this._kM = -1;

            this._time_attach = null;
            this._time_detach = null;
        }

        /// <summary>
        /// the constructor with parameter array
        /// </summary>
        /// <param name="_params">contains the parameter lists. 
        ///     the order is 0-ka1, 1-kd1, 2-ka2, 3-kd2, 4-ka3, 5-kd3, 6-conc, 7-Rmax, 8-R0_AB,
        ///     9-R0_AB_Star, 10-R0_B, 11-R0_B_Star
        ///     12-attachDuration, 13-detachDuration,14-deltaT</param>
        public ConformationSelection(double[] _params):base()
        {
            //check for the validity for the input
            if (_params.Count() < 14)
            {
                throw new Exception("the input parameter array is not valid. doesn't contain enough elements");
            }
            this._time_attach = null;
            this._time_detach = null;

            this._duration_attach = -1;
            this._duration_detach = -1;
            this._kM = -1;
            setParameters(_params);
            
            //this._fillTimeArrays();
        }

        /// <summary>
        /// here in this derived class, we also need to initialize the intermediate arrays,
        /// AB and AB_Start
        /// </summary>
        protected override void _fillTimeArrays()
        {
            base._fillTimeArrays();
            
            _R_AB_Att = new double[_ru_attach.Count];
            _R_AB_Star_Att = new double[_ru_attach.Count];

            _R_B_Att = new double[_ru_attach.Count];
            _R_B_Star_Att = new double[_ru_attach.Count];

            _R_AB_Det = new double[_ru_detach.Count];
            _R_AB_Star_Det = new double[_ru_detach.Count];

            _R_B_Det = new double[_ru_detach.Count];
            _R_B_Star_Det = new double[_ru_detach.Count];
            
            _R_AB_Att[0] = 0 ;
            _R_AB_Star_Att[0] = 0;

            _R_B_Att[0] = 0;//we set these two to be zero at this moment, but will change at the time we doing the simulation later
            _R_B_Star_Att[0] = 0;

            //these below things are set temporarily, will change later
            _R_AB_Det[0] = this.SSPR_r0;
            _R_AB_Star_Det[0] = this.SSPR_r0;
            _R_B_Det[0] = this.SSPR_r0;
            _R_B_Star_Det[0] = this.SSPR_r0;
        }

        /// <summary>
        /// now set the parameters!!! for duration, time, deltaT as well. for deltaT, it could have a default value. so we don't
        /// have to set in here, simply put -1 or omit it. Otherwise, we have to use the input.
        /// 
        /// there are two different way to do the calculation depending on R0.
        /// if we want to do the simulation, do we would not worry about R0, 
        /// let the system to take values from end of the association phase
        /// otherwise we will need to fit it as parameters. if we don't fit it
        /// we will put in negative values to indicate the system to take values
        /// from association phase end. 
        /// Normally, the R0s (such as R0_AB(*), R0_B(*)) is carried over from attaching phase.
        /// But for the real experiment, there are phase shift, we could get different value by the offset model.
        /// 
        /// it has the order, the order is 0-ka1, 1-kd1, 2-ka2, 3-kd2, 4-ka3, 5-kd3, 6-conc, 7-Rmax, 8-R0_AB,
        ///     9-R0_AB_Star, 10-R0_B, 11-R0_B_Star
        ///     12-attachDuration, 13-detachDuration,14-deltaT
        ///     So far DON"T have kM in there now.
        ///     here we also allow variable length of params, but only
        ///     allow default value of deltaT, duration of association or detach.
        ///     
        /// </summary>
        /// <param name="_params">it has the order, the order is 0-ka1, 1-kd1, 2-ka2, 3-kd2, 4-ka3, 5-kd3, 6-conc, 7-Rmax, 8-R0_AB,
        ///     9-R0_AB_Star, 10-R0_B, 11-R0_B_Star
        ///     12-attachDuration, 13-detachDuration,14-deltaT
        ///     So far DON"T have kM in there now.
        ///     here we also allow variable length of params, but only
        ///     allow default value of deltaT, duration of association or detach.
        ///     
        /// 
        ///     </param>
        public override void setParameters(double[] _params)
        {
            //here we allow variable length
            if (_params.Count() < 8) //8 is the minimum number of parameters to be set. since we are doing non-offset model. we can do this.
            {
                throw new Exception("the input parameter array is not valid. doesn't contain enough elements");
            }
            this._ka = _params[0];
            this._kd = _params[1];
            this._ka2 = _params[2];
            this._kd2 = _params[3];
            this._ka3 = _params[4];
            this._kd3 = _params[5];
            this._conc = _params[6];
            this._Rmax = _params[7];

            //now check for the extra parameters
            this.SSPR_r0 = _params[6];//by this we don't allow to remeber last round value
            if (_params.Count() >= 9 && _params[8] > 0)
            {
                //R0_AB = true;
                this.SSPR_r0 = _params[8];
            }
            
            this._R0_AB_Star = -1;
            if (_params.Count() >= 10 && _params[9] > 0)
            {
                this._R0_AB_Star = _params[9];
            }
            this._R0_B = -1;
            if (_params.Count() >= 11 && _params[10] > 0)
            {
                this._R0_B = _params[10];
            }
            this._R0_B_Star = -1;
            if (_params.Count() >= 12 && _params[11] > 0)
            {
                this._R0_B_Star = _params[11];
            }

            bool updateTimeArrays = false;//it could remeber last round value, that is the deal.
            if (_params.Count() >= 13 && _params[12] > 0)
            {
                updateTimeArrays = true;
                this._duration_attach = _params[12];
            }
            if (_params.Count() >= 14 && _params[13] > 0)
            {
                updateTimeArrays = true;
                this._duration_detach = _params[13];
            }
            if (_params.Count() >= 15 && _params[14] > 0)
            {
                updateTimeArrays = true;
                this._deltaT = _params[14];
            }

            if (updateTimeArrays)
            {
                _fillTimeArrays();
            }
        }
        /// <summary>
        /// NOT implemented in this class
        /// override the original one in the base class, but leave this one UNAVAILABLE in this two state model;
        /// </summary>
        /// <param name="_ka"></param>
        /// <param name="_kd"></param>
        /// <param name="_kM"></param>
        /// <param name="_conc"></param>
        /// <param name="_Rmax"></param>
        /// <param name="_r0"></param>
        public override void setParameters(double _ka, double _kd, double _kM, double _conc, double _Rmax, double _r0
                                                  )
        {
            throw new NotImplementedException("don't call this one, please");
        }

        /// <summary>
        /// a function used to add new values to the Time array
        /// </summary>
        /// <param name="_tVals">array of values to be added</param>
        public override void addValueSToTimeArr_attach(List<double> _tVals)
        {/*
            for (int i = 0; i < _tVals.Count; i++)
            {
                addValueToTimeArr_attach(_tVals[i]);
                if (i % 100 == 0)
                {
                    Console.WriteLine("doing i=" + i + " round of adding..............");
                }
            }
          
            if (_time_attach == null)
            {
                throw new Exception("unitialized time array");
            }
            _time_attach = MergeArrays(_time_attach, _tVals);*/
            //check the wether need to enlarge the AB and AB_star arrays, since they are new.
            base.addValueSToTimeArr_attach(_tVals);
            if (this._R_AB_Att.Count() != _time_attach.Count())
            {
                this._R_AB_Att = new double[_time_attach.Count()];
            }
            if (this._R_AB_Star_Att.Count() != _time_attach.Count())
            {
                this._R_AB_Star_Att = new double[_time_attach.Count()];
            }

            if (this._R_B_Att.Count() != _time_attach.Count())
            {
                this._R_B_Att = new double[_time_attach.Count()];
            }
            if (this._R_B_Star_Att.Count() != _time_attach.Count())
            {
                this._R_B_Star_Att = new double[_time_attach.Count()];
            }
        }



        public override void addValueSToTimeArr_detach(List<double> _tVals)
        {
            /*            if (_time_detach == null)
                        {
                            throw new Exception("unitialized time array");
                        }
            
                        for (int i = 0; i < _tVals.Count; i++)
                        {
                            addValueToTimeArr_detach(_tVals[i]);
                        }* /

                        _time_detach = MergeArrays(_time_detach, _tVals);
                        //check the wether need to enlarge the ru_attach
                        for (int i = _ru_detach.Count; i < _time_detach.Count(); i++)
                        {
                            _ru_detach.Add(0);
                        }
                        */
            base.addValueSToTimeArr_detach(_tVals);
            if (this._R_AB_Det.Count() != _time_detach.Count())
            {
                this._R_AB_Det = new double[_time_detach.Count()];
            }
            if (this._R_AB_Star_Det.Count() != _time_detach.Count())
            {
                this._R_AB_Star_Det = new double[_time_detach.Count()];
            }

            if (this._R_B_Det.Count() != _time_detach.Count())
            {
                this._R_B_Det = new double[_time_detach.Count()];
            }
            if (this._R_B_Star_Det.Count() != _time_detach.Count())
            {
                this._R_B_Star_Det = new double[_time_detach.Count()];
            }
        }

        /// <summary>
        /// this is euler sheme.
        /// </summary>
        public override void run_Attach()//this is the Euler
        {
            //check for the validity of input
            if (_ka < 0 || _kd < 0 || _ka2 < 0 || _kd2 < 0 || _ka3<0||_kd3<0||_conc < 0 || _Rmax < 0)
            {
                throw new Exception("unitialized parameters");
            }
            if (this._time_attach == null || this._time_detach == null)
            {
                throw new Exception("unitialized time arrays");
            }

            //start building the numerical integration
             //_ru.Add(0);the _ru_attach has been initialized and added with all zeros in the base class.
            _R_AB_Att[0] = 0 ;
            _R_AB_Star_Att[0] = 0;
            _ru_attach[0] = _R_AB_Star_Att[0] + _R_AB_Att[0];
            //in this conformation selection model, we need to determine R_B and R_B* at the beginning.
            //determine this by the ka3 and kd3 with a steady state
            this._R_B_Star_Att[0] = this._Rmax / (_kd3 / _ka3 + 1);
            _R_B_Att[0] = _Rmax - this._R_B_Star_Att[0];

            //now star doing the differential equation
            double dR_AB_Star_per_dt;
            double dR_AB_per_dt;
            double dR_B_Star_per_dt, dR_B_per_dt;
            double dt;
            for( int i=0; ;i++)
            {
                double R_B_i = _Rmax - _R_B_Star_Att[i] - _R_AB_Att[i] - _R_AB_Star_Att[i];
	            dR_AB_per_dt=_ka*_conc*R_B_i
                                        - _kd * _R_AB_Star_Att[i];// +_kd2 * _R_AB_Att[i] - _ka2 * _R_AB_Star_Att[i];
                dR_AB_Star_per_dt = _ka2 * _conc*_R_B_Star_Att[i] - _kd2 * _R_AB_Star_Att[i];

                dR_B_per_dt = -1 * dR_AB_per_dt - _ka3 * R_B_i + _kd3 * _R_B_Star_Att[i];
                dR_B_Star_per_dt = -1 * dR_AB_Star_per_dt + _ka3 * R_B_i - _kd3 * _R_B_Star_Att[i];
	            
	            if(i>=_ru_attach.Count -1 )
	                {
                        break;
		                
	                }
                dt = this._time_attach[i + 1] - this._time_attach[i];
                _R_AB_Att[i + 1] = dR_AB_per_dt * dt + _R_AB_Att[i];
                _R_AB_Star_Att[i + 1] = dR_AB_Star_per_dt * dt + _R_AB_Star_Att[i];
                _R_B_Star_Att[i + 1] = dR_B_Star_per_dt * dt + _R_B_Star_Att[i];
                _R_B_Att[i + 1] = _Rmax - _R_B_Star_Att[i + 1] - _R_AB_Att[i + 1] - _R_AB_Star_Att[i + 1];
                _ru_attach[i + 1] = _R_AB_Att[i+1]+_R_AB_Star_Att[i+1];

            }
        
        }
        public  void run_Attach_RK()//this is the Runge-Kutta
        {
            //check for the validity of input
            if (_ka < 0 || _kd < 0 || _ka2 < 0 || _kd2 < 0 ||_ka3<0||_kd3<0|| _conc < 0 || _Rmax < 0)
            {
                throw new Exception("unitialized parameters");
            }
            if (this._time_attach == null || this._time_detach == null)
            {
                throw new Exception("unitialized time arrays");
            }

            //start building the numerical integration
            //_ru.Add(0);the _ru_attach has been initialized and added with all zeros in the base class.
            _R_AB_Att[0] = 0;
            _R_AB_Star_Att[0] = 0;
            _ru_attach[0] = _R_AB_Star_Att[0] + _R_AB_Att[0];
            List<List<double>> RABs = new List<List<double>>();
            //initialize the holder array for RABs
            for (int i = 0; i < _R_AB_Att.Count(); i++)
            {
                RABs.Add(new List<double>{_R_AB_Star_Att[i],_R_AB_Att[i], _R_B_Star_Att[i], _R_B_Att[i]});
            }
           
            //RUN runge-kutta 4th order.
            RungeKuttaMethod.RungeKutta.Solution(this.Derivative_Attach, this._time_attach, ref RABs);

            //now let's unpack the output
            for (int i = 1; i < RABs.Count; i++)
            {
                this._R_AB_Att[i] = RABs[i][1];
                this._R_AB_Star_Att[i] = RABs[i][0];
                this._R_B_Star_Att[i] = RABs[i][2];
                this._R_B_Att[i] = RABs[i][3];
                this._ru_attach[i] = _R_AB_Att[i] + _R_AB_Star_Att[i];
            }
            //done
        }


        /// <summary>
        /// this is the derivative function or the function for differential equation 
        /// to do the RK or regular 
        /// </summary>
        /// <param name="_t"></param>
        /// <param name="_RABs"></param>
        /// <returns></returns>
        protected List<double> Derivative_Attach(double _t, List<double> _RABs)
        {
            List<double> nextRABs = new List<double>();
            
            //in this _RABs, 0 is for R_AB_Star and 1 is for R_AB, 2 for R_B_Star, 3 for R_B
            //double R_B_i = _Rmax - _RABs[0] - _RABs[1] - _RABs[2];
            double dR_AB_Star_per_dt = _ka2 * _conc * _RABs[2]
                                        - _kd2 * _RABs[0];// +_kd2 * _RABs[1] - _ka2 * _RABs[0];
            double dR_AB_per_dt = _ka *_conc* _RABs[3] - _kd * _RABs[1];

            double dR_B_Star_per_dt = -1*dR_AB_Star_per_dt+_ka3 * _RABs[3] - _kd3 * _RABs[2];

            nextRABs.Add(dR_AB_Star_per_dt);
            nextRABs.Add(dR_AB_per_dt);
            nextRABs.Add(dR_B_Star_per_dt);
            nextRABs.Add(0-dR_B_Star_per_dt-dR_AB_per_dt-dR_AB_Star_per_dt); //the fact is that dB+dB*+dAB+dAB*=0; all changes sum to zero.
            return nextRABs;
        }

        public override void run_Detach() //Euler scheme
        {
            //check for the validity of input
            if (_ka < 0 || _kd < 0 || _ka2 < 0 || _kd2 < 0 || _ka3<0||_kd3<0||_conc < 0 || _Rmax < 0)
            {
                throw new Exception("unitialized parameters");
            }
            if (this._time_attach == null || this._time_detach == null)
            {
                throw new Exception("unitialized time arrays");
            }
            
            //in this one, it is abit tricky, we need to figure out what to do with R0
            if (SSPR_r0 < 0 && _R0_AB_Star < 0)
            {
                this._R_AB_Det[0] = _R_AB_Att[_R_AB_Att.Count() - 1];
                this._R_AB_Star_Det[0] = _R_AB_Star_Att[_R_AB_Star_Att.Count() - 1];
            }
            else
            {
                this._R_AB_Det[0] = this.SSPR_r0 ;
                this._R_AB_Star_Det[0] = this._R0_AB_Star;
            }
            
            //for R_B and R_B*
            if (this._R0_B < 0 && _R0_B_Star < 0)
            {
                this._R_B_Det[0] = _R_B_Att[_R_B_Att.Count() - 1];
                this._R_B_Star_Det[0] = _R_B_Star_Att[_R_B_Star_Att.Count() - 1];
            }
            else
            {
                this._R_B_Det[0] = this._R0_B ;
                this._R_B_Star_Det[0] = this._R0_B_Star;
            }
            
            _ru_detach[0] = _R_AB_Det[0] + _R_AB_Star_Det[0];

            //now doing the calculation
            double dR_AB_Star_per_dt;
            double dR_AB_per_dt;
            double dR_B_Star_per_dt, dR_B_per_dt;
            double dt;
            //double conc0 = 0;
            for (int i = 0; ; i++)
            {
                double R_B_i = _Rmax - _R_B_Star_Det[i] - _R_AB_Det[i] - _R_AB_Star_Det[i];
                dR_AB_per_dt = 0//_ka * _conc * R_B_i
                                        - _kd * _R_AB_Star_Det[i];// +_kd2 * _R_AB_Att[i] - _ka2 * _R_AB_Star_Att[i];
                dR_AB_Star_per_dt = 0//_ka2 * _conc * _R_B_Star_Att[i] 
                        - _kd2 * _R_AB_Star_Det[i];

                dR_B_per_dt = -1 * dR_AB_per_dt - _ka3 * R_B_i + _kd3 * _R_B_Star_Det[i];
                dR_B_Star_per_dt = -1 * dR_AB_Star_per_dt + _ka3 * R_B_i - _kd3 * _R_B_Star_Det[i];

                if (i >= _ru_detach.Count - 1)
                {
                    break;

                }
                dt = this._time_detach[i + 1] - this._time_detach[i];
                _R_AB_Det[i + 1] = dR_AB_per_dt * dt + _R_AB_Det[i];
                _R_AB_Star_Det[i + 1] = dR_AB_Star_per_dt * dt + _R_AB_Star_Det[i];
                _R_B_Star_Det[i + 1] = dR_B_Star_per_dt * dt + _R_B_Star_Det[i];
                _R_B_Det[i + 1] = _Rmax - _R_B_Star_Det[i + 1] - _R_AB_Det[i + 1] - _R_AB_Star_Det[i + 1];
                _ru_detach[i + 1] = _R_AB_Det[i + 1] + _R_AB_Star_Det[i + 1];
                
            }
            //done!!
        }

        public void run_Detach_RK() //Runge-Kutta scheme
        {
            //check for the validity of input
            if (_ka < 0 || _kd < 0 || _ka2 < 0 || _kd2 < 0 ||_ka3<0||_kd3<0|| _conc < 0 || _Rmax < 0)
            {
                throw new Exception("unitialized parameters");
            }
            if (this._time_attach == null || this._time_detach == null)
            {
                throw new Exception("unitialized time arrays");
            }

            //in this one, it is abit tricky, we need to figure out what to do with R0
            if (SSPR_r0 < 0 && _R0_AB_Star < 0)
            {
                this._R_AB_Det[0] = _R_AB_Att[_R_AB_Att.Count() - 1];
                this._R_AB_Star_Det[0] = _R_AB_Star_Att[_R_AB_Star_Att.Count() - 1];
            }
            else
            {
                this._R_AB_Det[0] = this.SSPR_r0;
                this._R_AB_Star_Det[0] = this._R0_AB_Star;
            }

            //for R_B and R_B*
            if (this._R0_B < 0 && _R0_B_Star < 0)
            {
                this._R_B_Det[0] = _R_B_Att[_R_B_Att.Count() - 1];
                this._R_B_Star_Det[0] = _R_B_Star_Att[_R_B_Star_Att.Count() - 1];
            }
            else
            {
                this._R_B_Det[0] = this._R0_B;
                this._R_B_Star_Det[0] = this._R0_B_Star;
            }
            _ru_detach[0] = _R_AB_Det[0] + _R_AB_Star_Det[0];

            List<List<double>> RABs = new List<List<double>>();
            for (int i = 0; i < _R_AB_Det.Count(); i++)
            {
                RABs.Add(new List<double> { _R_AB_Star_Det[i], _R_AB_Det[i], _R_B_Star_Det[i], _R_B_Det[i] });
            }
            //doing the Runge-kutta
            RungeKuttaMethod.RungeKutta.Solution(this.Derivative_Detach, this._time_detach, ref RABs);

            //now let's unpack the output
            for (int i = 1; i < RABs.Count; i++)
            {
                this._R_AB_Det[i] = RABs[i][1];
                this._R_AB_Star_Det[i] = RABs[i][0];
                this._R_B_Star_Det[i] = RABs[i][2];
                this._R_B_Det[i] = RABs[i][3];
                this._ru_detach[i] = _R_AB_Det[i] + _R_AB_Star_Det[i];
            }
            //done!!
        }

        protected List<double> Derivative_Detach(double _t, List<double> _RABs)
        {
            List<double> nextRABs = new List<double>();

            //in this _RABs, 0 is for R_AB_Star and 1 is for R_AB
            //double dR_AB_Star_per_dt = 0 //_ka * _conc * (_Rmax - _RABs[1] - _RABs[0])
                                        //- _kd * _RABs[0] + _kd2 * _RABs[1] - _ka2 * _RABs[0];
            //double dR_AB_per_dt = _ka2 * _RABs[0] - _kd2 * _RABs[1];
            double dR_AB_Star_per_dt = 0//_ka2 * _conc * _RABs[2]
                                        - _kd2 * _RABs[0];// +_kd2 * _RABs[1] - _ka2 * _RABs[0];
            double dR_AB_per_dt = 0//_ka *_conc* _RABs[3] 
                                    - _kd * _RABs[1];

            double dR_B_Star_per_dt=-1*dR_AB_Star_per_dt +_ka3*_RABs[3]-_kd3*_RABs[2];

            nextRABs.Add(dR_AB_Star_per_dt);
            nextRABs.Add(dR_AB_per_dt);
            nextRABs.Add(dR_B_Star_per_dt);
            nextRABs.Add(0 - dR_B_Star_per_dt - dR_AB_Star_per_dt - dR_AB_per_dt);
            return nextRABs;
        }

 

        public double[] RU_AB_Attach
        {
            get { return _R_AB_Att; }
        }
        public double[] RU_AB_Star_Attach
        {
            get { return _R_AB_Star_Att; }
        }
        public double[] RU_AB_Detach
        {
            get { return _R_AB_Det; }
        }
        public double[] RU_AB_Star_Detach
        {
            get { return _R_AB_Star_Det; }
        }

        public double[] RU_B_Attach
        {
            get { return _R_B_Att; } 
        }
        public double[] RU_B_Star_Attach
        {
            get { return this._R_B_Star_Att ;}
        }
        public double[] RU_B_Detach
        {
            get { return _R_B_Det; }
        }
        public double[] RU_B_Star_Detach
        {
            get { return this._R_B_Star_Det; }
        }
        
    }//end of class
}
