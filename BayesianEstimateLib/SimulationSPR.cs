﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


using Meta.Numerics;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;
using AccessoryLib;

namespace BayesianEstimateLib
{
    
    public abstract class SimulationSPR
    {
        const double EPSILON = 1E-10;
        /* don't allow empty constructor in order to keep the integrity of data
        /// <summary>
        /// empty constructer
        /// </summary>
        public simulationSPR()
        {
            //empty constructer
        }
         * */
        /// <summary>
        /// the constructer with necessary parameters
        /// </summary>
        /// <param name="_ka">on rate constant, unit 1/m*1/s, 1E4 - 1E6</param>
        /// <param name="_kd">off rate constant, unit -/s, 1E-3 - 1E-6 (??)</param>
        /// <param name="conc">the concentration of analytes in the flow buffer, unit M (molar/l)</param>
        /// <param name="Rmax">the maximum Response Unit, unit arbitrary unit</param>
        /// <param name="duration">the association duration, unit second</param>
        /// <param name="kM">the diffusion rate, unit m/s, about 1E-5 in the BIAcore flow cell.</param>
        /// <param name="deltaT">time step to run numerical integration, the smaller the better, 1/100s or smaller</param>
        public SimulationSPR(double _ka, double _kd, double _conc, double _Rmax, double _r0,
                                                double _duration_attach, double _duration_detach,  double _kM=1E6, double _deltaT=0.01 )
        {
            this._ka = _ka;
            this._kd = _kd;
            this._conc = _conc;
            this._Rmax = _Rmax;
            this._duration_attach = _duration_attach;
            this._duration_detach = _duration_detach;
            this._kM = _kM;
            this._deltaT = _deltaT;
            this.SSPR_r0 = _r0;
            
            _time_attach = new List<double>();
            _time_detach = new List<double>();
            _ru_attach = new List<double>();
            _ru_detach = new List<double>();
            for (int i = 0; i * _deltaT < _duration_attach; i++)
            {
                _time_attach.Add(i * _deltaT);
                _ru_attach.Add(0);
            }
            for (int i = 0; i * _deltaT < _duration_detach; i++)
            {
                _time_detach.Add(i * _deltaT);
                _ru_detach.Add(0);
            }
        }

        public void setParameters(double _ka, double _kd, double _kM, double _conc, double _Rmax, double _r0
                                                  )
        {
            this._ka = _ka;
            this._kd = _kd;
            this._kM = _kM;
            this._conc = _conc;
            this._Rmax = _Rmax;
            this.SSPR_r0 = _r0;
        }

        /// <summary>
        /// used to generate normal distributed noise to the RU 
        /// </summary>
        public void addNoise(double _sigma=0.05)
        {
            NormalDistribution normGenerateor = new NormalDistribution(0, _sigma);
            
            Random rng = new Random(AccessoryLib.AceessoryLib.SEED);

            
            for (int i = 0; i < _ru_attach.Count; i++)
            {
                _ru_attach[i] += normGenerateor.GetRandomValue(rng);
            }

            for (int i = 0; i < _ru_detach.Count; i++)
            {
                _ru_detach[i] += normGenerateor.GetRandomValue(rng);
            }
        }


        public abstract void run_Attach();
        
        public abstract void run_Detach();
        

        /// <summary>
        /// a function used to add new values to the Time array
        /// </summary>
        /// <param name="_tVals">array of values to be added</param>
        public void addValueSToTimeArr_attach(List<double> _tVals)
        {/*
            for (int i = 0; i < _tVals.Count; i++)
            {
                addValueToTimeArr_attach(_tVals[i]);
                if (i % 100 == 0)
                {
                    Console.WriteLine("doing i=" + i + " round of adding..............");
                }
            }
          */
            _time_attach = MergeArrays(_time_attach, _tVals);
            //check the wether need to enlarge the ru_attach
            for (int i = _ru_attach.Count; i < _time_attach.Count(); i++)
            {
                _ru_attach.Add(0);
            }

        }

        

        public void addValueSToTimeArr_detach(List<double> _tVals)
        {
            /*
            for (int i = 0; i < _tVals.Count; i++)
            {
                addValueToTimeArr_detach(_tVals[i]);
            }*/

            _time_detach = MergeArrays(_time_detach, _tVals);
            //check the wether need to enlarge the ru_attach
            for (int i = _ru_detach.Count; i < _time_detach.Count(); i++)
            {
                _ru_detach.Add(0);
            }

        }
        /*
        /// <summary>
        /// in this method, we assume the numerical solution has been runned and we just picked the appropriated values based on the input t values
        /// also we assume the values has been added before the simulations.
        /// </summary>
        /// <param name="_tVals"></param>
        /// <returns></returns>
        public List<double> GetRUValuesAttach(List<double> _tVals)
        {
            List<double> fittedVals = new List<double>(_tVals.Count);
            for (int i = 0; i < _tVals.Count; i++)
            {
                //need to find the specific t value in the numerical solution
                for (int j = 0; j < this._time_attach.Count; j++)
                {
                    if ((this._time_attach[j] - _tVals[i]) < EPSILON && (this._time_attach[j] - _tVals[i]) > -1 * EPSILON)
                        fittedVals.Add(this._ru_attach[j]);
                }
            }
            if (fittedVals.Count != _tVals.Count)
                throw new System.Exception("not equal size array when finding values. something wrong");
            return fittedVals;
        }*/

        /// <summary>
        /// in this method, we assume the numerical solution has been runned and we just picked the appropriated values based on the input t values
        /// also we assume the values has been added before the simulations. Now this is implemented with MergeSort-like algorithm 
        /// </summary>
        /// <param name="_tVals"></param>
        /// <returns></returns>
        public List<double> GetRUValuesAttach(List<double> _tVals)
        {
            List<double> fittedVals = new List<double>(_tVals.Count);
            int indexAiming = 0, indexHolding = 0;//two index to loop through the two arrays to find the specific items
            while (true)
            {   
                //we need to check whether we are going ahead, 
                if (indexAiming >= _tVals.Count || indexHolding > this._time_attach.Count)
                    break;//for either one is finished,we are done
                
                if ((this._time_attach[indexHolding] - _tVals[indexAiming]) < EPSILON && (this._time_attach[indexHolding] - _tVals[indexAiming]) > -1 * EPSILON)
                {//find one matching
                    fittedVals.Add(this._ru_attach[indexHolding]);
                    indexAiming++;
                    indexHolding++;
                    //increment both index, we assume there is no tie in either of the array
                }
                else
                {
                    //there is no matching, so we only increment the holding array to look at the next one in holding array.
                    //keep the current aiming index the same
                    indexHolding++;
                }
                
            }

            //now we are out, so we need to check whether we are doing a good job,
            if (indexAiming != _tVals.Count)
            {
                //now something wrong, we need to throw exceptions
                throw new System.Exception("not equal size array when finding values. something wrong");
            }
            //if the indexHoling is not the end, if could be OK.
            /*
            for (int i = 0; i < _tVals.Count; i++)
            {
                //need to find the specific t value in the numerical solution
                for (int j = 0; j < this._time_attach.Count; j++)
                {
                    if ((this._time_attach[j] - _tVals[i]) < EPSILON && (this._time_attach[j] - _tVals[i]) > -1 * EPSILON)
                        fittedVals.Add(this._ru_attach[j]);
                }
            }
            if (fittedVals.Count != _tVals.Count)
                throw new System.Exception("not equal size array when finding values. something wrong");
             * */
            return fittedVals;
        }
        /*
        /// <summary>
        /// in this method, we assume the numerical solution has been runned and we just picked the appropriated values based on the input t values
        /// also we assume the values has been added before the simulations.
        /// </summary>
        /// <param name="_tVals"></param>
        /// <returns></returns>
        public List<double> GetRUValuesDetach(List<double> _tVals)
        {
            List<double> fittedVals = new List<double>(_tVals.Count);
            for (int i = 0; i < _tVals.Count; i++)
            {
                //need to find the specific t value in the numerical solution
                for (int j = 0; j < this._time_detach.Count; j++)
                {
                    if ((this._time_detach[j] - _tVals[i]) < EPSILON && (this._time_detach[j] - _tVals[i]) > -1 * EPSILON)
                        fittedVals.Add(this._ru_detach[j]);
                }
            }
            if (fittedVals.Count != _tVals.Count)
                throw new System.Exception("not equal size array when finding values. something wrong");
            return fittedVals;
        }*/
        /// <summary>
        /// in this method, we assume the numerical solution has been runned and we just picked the appropriated values based on the input t values
        /// also we assume the values has been added before the simulations. using mergeSort algorithm
        /// </summary>
        /// <param name="_tVals"></param>
        /// <returns></returns>
        public List<double> GetRUValuesDetach(List<double> _tVals)
        {
            List<double> fittedVals = new List<double>(_tVals.Count);
            int indexAiming = 0, indexHolding = 0;//two index to loop through the two arrays to find the specific items
            while (true)
            {
                //we need to check whether we are going ahead, 
                if (indexAiming >= _tVals.Count || indexHolding > this._time_detach.Count)
                    break;//for either one is finished,we are done

                if ((this._time_detach[indexHolding] - _tVals[indexAiming]) < EPSILON && (this._time_detach[indexHolding] - _tVals[indexAiming]) > -1 * EPSILON)
                {//find one matching
                    fittedVals.Add(this._ru_detach[indexHolding]);
                    indexAiming++;
                    indexHolding++;
                    //increment both index, we assume there is no tie in either of the array
                }
                else
                {
                    //there is no matching, so we only increment the holding array to look at the next one in holding array.
                    //keep the current aiming index the same
                    indexHolding++;
                }

            }

            //now we are out, so we need to check whether we are doing a good job,
            if (indexAiming != _tVals.Count)
            {
                //now something wrong, we need to throw exceptions
                throw new System.Exception("not equal size array when finding values. something wrong");
            }
            return fittedVals;
        }

        /// <summary>
        /// the merge method to combine the two array into one, contains the union of the two array. this method assumes the two arrays are sorted in ascending order
        /// </summary>
        /// <param name="_a1">input array one</param>
        /// <param name="_a2">input array two</param>
        /// <returns>the return union array</returns>
        public static List<double> MergeArrays(List<double> _a1, 
            List<double> _a2)
        {
            double epsilon=1E-8;
            List<double> ret=new List<double>();
            int a1_index = 0, a2_index = 0;
            while (true)
            {
                if ((_a1[a1_index] - _a2[a2_index]) < epsilon && (_a1[a1_index] - _a2[a2_index]) > -1 * epsilon)
                {
                    //two values equal
                    //then we get any of it and also move forward for both
                    ret.Add(_a1[a1_index]);
                    a1_index++;
                    a2_index++;
                }
                else //not equal, we need to add the smaller one and move forward.
                {
                    if (_a1[a1_index] > _a2[a2_index])
                    {
                        ret.Add(_a2[a2_index]);
                        a2_index++;
                    }
                    else
                    {
                        ret.Add(_a1[a1_index]);
                        a1_index++;
                    }
                }

                //now we need to check wether a1 or a2 has been finished.
                if (a1_index >= _a1.Count || a2_index >= _a2.Count)
                {
                    break;
                }
            }//end of while

            //now we need to copy over the left over for one array.
            for (; a1_index < _a1.Count; a1_index++)
            {
                ret.Add(_a1[a1_index]);
            }
            for(; a2_index < _a2.Count;a2_index++ )
            {
                ret.Add(_a2[a2_index]);
            }

            return ret;
        }
        //members
        //check the definition for these paramter in the constructor definition above
        protected double _ka;
        protected double _kd;
        protected double _conc;
        protected double _Rmax;
        protected double _duration_attach;
        protected double _duration_detach;
        protected double _kM;
        protected double _deltaT;
        protected double SSPR_r0;
        
        protected List<double> _time_attach;
        protected List<double> _time_detach;
        protected List<double> _ru_attach;
        protected List<double> _ru_detach;

        //properties
        /*public double ka
        {
            set { _ka = value; }
            get { return _ka; }
        }

        public double kd
        {
            set { _kd = value; }
            get { return _kd; }
        }

        public double Conc
        {
            set { _conc = value; }
            get { return _conc; }
        }

        public double Rmax
        {
            set { _Rmax = value; }
            get { return _Rmax; }
        }
        public double kM
        {
            set { _kM = value; }
            get { return _kM; }
        }
        */
        public double Duration_Attach
        {
            //set { _duration_attach = value; }
            get { return _duration_attach; }
        }

        public double _Duration_Detach
        {
            set { _duration_detach = value; }
            get { return _duration_detach; }
        }

        public double DeltaT
        {
            set { _deltaT = value; }
            get { return _deltaT; }
        }

        public List<double> Time_Detach
        {
            set { _time_detach = value; }
            get { return _time_detach; }
        }
        public List<double> Time_Attach
        {
            set { _time_attach = value; }
            get { return _time_attach; }
        }
        public List<double> RU_Detach
        {
            //set { _ru = value; }
            get { return _ru_detach; }
        }
        public double  R0
        {
            set { SSPR_r0 = value; }
            get { return SSPR_r0 ; }
        }
        public List<double> RU_Attach
        {
            //set { _ru = value; }
            get { return _ru_attach; }
        }
    }//end of class
}
