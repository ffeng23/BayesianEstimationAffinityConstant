using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
//using BayesianEstimateLib;
namespace BayesianEstimateLib
{
    /*
    /// <summary>
    /// here this is a delegate, with _x as time array and parameters. return a series of data
    /// like the two state or any other model to calculate the SPR sensorgram
    /// </summary>
    /// <param name="_params"></param>
    /// <param name="_x"></param>
    /// <returns></returns>
    public delegate double[] ModelFuncDelegate(double[] _params,double _conc, double[] _ts);
    */
    /// <summary>
    /// this class is the data structur to hold the SPR sensorgram data
    /// data don't know the model and they don't have to know the model
    /// </summary>
    public class SPR_Data
    {
        //this is the constructor
        public SPR_Data()
        {
            //empty one
            this._FlagConcRead = false;
            this._NumOfSeries = -1;
            this._SPRDataSeries = null;
            /*
            this._NumOfCommonParams = _numOfCommonParam;
            this._NumOfIndividualParams = _numOfIndividualParam;
            this._Common_Parameters = new double[this._NumOfCommonParams];
            this._ModelFuncDel = _funcDel;
            */
        }
        
        /// <summary>
        /// at this point, we assume the association time interval and dissociation time interval are reasonal 
        /// no overlaping and association is before the dissociation time and start is smaller than end;
        /// </summary>
        /// <param name="_FileName"></param>
        /// <param name="ass_start"></param>
        /// <param name="ass_end"></param>
        /// <param name="dis_start"></param>
        /// <param name="dis_end"></param>
        public void ReadDataFromDict(Dictionary<int, List<double>> dict, double ass_start , double ass_end, 
            double dis_start, double dis_end )
        {
            //this._FileName = _FileName;
            //Dictionary<int, List<double> > dict=BayesianEstimateLib.DataIO.ReadDataTable(_FileName);
            
            //first we need check for consistency.
            if (dict.Count/2*2!=dict.Count||dict.Count/2 != this._Concs.Count())
            {
                throw new Exception("the number of concentrations and the number of data series are not consistent. please check");
            }

            //now we are good, just read the data into the variables.
            this._SPRDataSeries = new SPR_DataSerie[_Concs.Count()];
            for (int i=0;i<_Concs.Count();i++)
            {
                _SPRDataSeries[i] = new SPR_DataSerie(/*this._ModelFuncDel*/);
                _SPRDataSeries[i].SetConc(this._Concs[i]);
                List<double> ass_time=new List<double>();
                List<double> ass_ru=new List<double>();
                List<double> dis_time=new List<double>();
                List<double> dis_ru=new List<double>();
                
                for(int j=0;j<dict[i*2].Count();j++)
                {
                    if (Double.IsNaN(dict[i*2][j]) || Double.IsNaN(dict[i*2+1][j]))
                    {
                        continue;
                    }
                    //check to see which one is go witch.
                    if (dict[i*2][j] >= ass_start && dict[i*2][j] <= ass_end)
                    {
                        //association phase data
                        ass_time.Add(dict[i*2][j]);
                        ass_ru.Add(dict[i*2 + 1][j]);
                    }
                    if (dict[i*2][j] >= dis_start && dict[i*2][j] <= dis_end)
                    {
                        //dissociation phase data
                        dis_time.Add(dict[i*2][j]-dis_start);
                        dis_ru.Add(dict[i*2 + 1][j]);
                    }
                }
                if (dis_time.Count() == 0 | ass_time.Count() == 0)
                {
                    throw new Exception("no data points read in for association/dissociation phase. please check");
                }
                //ass_time.Sort();
                
                this._SPRDataSeries[i].SetNumberOfAssociationPoints(ass_time.Count());
                this._SPRDataSeries[i].SetNumberOfDissociationPoints(dis_time.Count());
                this._SPRDataSeries[i].SetTimeOfAssociation(ass_time);
                this._SPRDataSeries[i].SetTimeOfDissociation(dis_time);
                this._SPRDataSeries[i].SetRUAssociation(ass_ru);
                this._SPRDataSeries[i].SetRUDissociation(dis_ru);

                //now we want to do the sorting to put the thing in order.
                int[] idx;// = new double[ass_time.Count()];
                _FillIndexArray(out idx, ass_time.Count());
                Array.Sort(this._SPRDataSeries[i].AssociationTimes , idx);
                //now we get the index sorted. so rearrange
                this._SPRDataSeries[i].SetRUAssociation(new List<double>(_RearrangeArray(this._SPRDataSeries[i].AssociationRUs, idx)));

                //now we want to do the dissiciation part
                _FillIndexArray(out idx, ass_time.Count());
                Array.Sort(this._SPRDataSeries[i].DissociationTimes, idx);
                //now we get the index sorted. so rearrange
                this._SPRDataSeries[i].SetRUDissociation(new List<double>(_RearrangeArray(this._SPRDataSeries[i].DissociationRUs, idx)));
                //done, but need to test this part!!!
            }

            //done
        }
        private void _FillIndexArray(out int[] _idx, int _count)
        {
            _idx = new int[_count];
            for (int i = 0; i < _count; i++)
            {
                _idx[i] = i;
            }
        }

        /// <summary>
        /// rearrange the elements in the array arrording to index array _idx
        /// </summary>
        /// <param name="_array">data array to be rearranged, this is also an output</param>
        /// <param name="_idx">the index array that holding the correct index for the output array</param>
        private double[] _RearrangeArray(double[] _array, int[] _idx)
        {
            //now first we need to reserve some space to hold the temporary values
            //is there a better algorithm that doesn't need any temporary space
            double[] temp = new double[_array.Count()];

            //now doing the rearranging, and only use the temp space where necessary
            //averagely, it used half of the space for rearranging
            for (int i = 0; i < _array.Count(); i++)
            {
                    temp[i] = _array[_idx[i]];
            }

            //done!!!
            return temp;
            
        }
        //
        public string ReadConcentrationFile(string _FileName)
        {
            string ConcStr="";
            //first open up the file
            string line;
            //int nRow = 0;
            //Dictionary<int, List<Double>> output = new Dictionary<int, List<double>>();
            FileInfo src = new FileInfo(_FileName);
            if (!src.Exists)
            {
                throw new System.Exception("file doesn't exist");
                // Environment.Exit(-1);

            }
            TextReader rin = src.OpenText();
            line = rin.ReadLine();
            
            while (line != null)
            {
                ConcStr = ConcStr+ line+"\r\n";
                line = rin.ReadLine();
            }
            rin.Close();
            _FlagConcRead = false;
            ReadConcentrationTextBox(ConcStr);
            return ConcStr;
        }
        public void ReadConcentrationTextBox(string _cstr)
        {
            if(this._FlagConcRead)
                return ;
            string[] temp = _cstr.Split(new char[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
            this._Concs = new double[temp.Count()];

            for (int i = 0; i < temp.Count(); i++)
            {
                _Concs[i] = Convert.ToDouble(temp[i]);
            }
            _FlagConcRead = true;
            this._NumOfSeries = temp.Count();
        }
        public void ResetFlagConcRead()
        {
            this._FlagConcRead = false;
        }

        /*
        /// <summary>
        /// calculate
        /// </summary>
        /// <param name="_parameters"></param>
        /// <returns></returns>
        public double LogLikelihood(double[] _parameters, ModelFuncDelegate _ModelFuc)
        {
            double ll = 0;
            for (int i = 0; i < NSeries; i++)
            {
                // repack parameter vector for each dilution series
                double[] parametersShort = new double[this.C_NumOfCommonParams + this.C_NumOfIndividualParams];
                for (int j = 0; j < this.C_NumOfIndividualParams; j++)
                {
                    parametersShort[j] = _parameters[i * this.C_NumOfIndividualParams + j];
                }
                for (int k = 0; k < this.C_NumOfCommonParams; k++)
                {
                    parametersShort[this.C_NumOfIndividualParams + k] = _parameters[NSeries * this.C_NumOfIndividualParams + k];
                }
                //                parametersShort[2] = parameters[parameters.Length - 1]; // beta
                ll += DS[i].LogLikelihood(parametersShort);
            }
            return ll;
        }
        */

        public int NSeries
        {
            get { return _NumOfSeries; }
        }
        public SPR_DataSerie[] SPRDataSeries
        {
            get { return _SPRDataSeries; }
        }
        private int _NumOfSeries;//this is num of data series like different conc
        private double[] _Concs; //it is possible that each 
        private SPR_DataSerie[] _SPRDataSeries;
        private bool _FlagConcRead;//used to indicated whether the concentration has been read in
            //if not we need to read the text field. we will never need to read the file  again.
            //the file can be called to read a new by the user, but not by the program
            //the user might have type in/modify the text fields, so we might need to read again the
            //text fields.

        private string _FileName;
        public string FileName
        {
            get { return this._FileName; }
        }
        /*protected int _NumOfIndividualParams;
        protected int _NumOfCommonParams;
        public int NumOfCommonParameters
        {
            get { return this._NumOfCommonParams; }
        }
        public int NumOfInidividualParameters
        {
            get { return this._NumOfIndividualParams; }
        }
        protected double[] _Common_Parameters;
        public double[] Common_Parameters
        {
            get { return this._Common_Parameters; }
            set { this._Common_Parameters = value; }
        }

        protected ModelFuncDelegate _ModelFuncDel;*/
    }//end of class
    public class SPR_DataSerie
    {
        public SPR_DataSerie(/*double[] _params_att, ModelFuncDelegate _mfuc_att, ModelFuncDelegate _mfuc_det*/)
        {
            this._NumOfDissociationPoints = -1;
            this._NumOfAssociationPoints = -1;
        }

        public void SetNumberOfAssociationPoints(int _N)
        {
                this._NumOfAssociationPoints=_N;
        }
        public void SetNumberOfDissociationPoints(int _N)
        {
            this._NumOfDissociationPoints=_N;
        }

        public void SetTimeOfAssociation(List<double> _t)
        {
            if (this._NumOfAssociationPoints == -1)
            {
                this._NumOfAssociationPoints = _t.Count;
            }
            //we assume the input array is well prepared and we don't check for errors
            if (this._NumOfAssociationPoints !=-1&&_t.Count != this._NumOfAssociationPoints)
            {
                throw new Exception("array are not set correctly with not correct number of elements");
            }
            _TimeAssociation = _t.ToArray ();//it is not very important to be in order ???
        }
        public void SetTimeOfDissociation(List<double> _t)
        {
            if (this._NumOfDissociationPoints == -1)
            {
                this._NumOfDissociationPoints = _t.Count;
            }
            //we assume the input array is well prepared and we don't check for errors
            if (this._NumOfDissociationPoints != -1 && _t.Count != this._NumOfDissociationPoints)
            {
                throw new Exception("array are not set correctly with not correct number of elements");
            }
            _TimeDissociation = _t.ToArray();//it is not very important to be in order ???
        }
        public void SetRUAssociation(List<double> _r)
        {
            if (this._NumOfAssociationPoints == -1)
            {
                this._NumOfAssociationPoints = _r.Count;
            }
            //we assume the input array is well prepared and we don't check for errors
            if (this._NumOfAssociationPoints != -1 && _r.Count != this._NumOfAssociationPoints)
            {
                throw new Exception("array are not set correctly with not correct number of elements");
            }
            _RUAssociation = _r.ToArray();//it is not very important to be in order ???
        }
        public void SetRUDissociation(List<double> _r)
        {
            if (this._NumOfDissociationPoints == -1)
            {
                this._NumOfDissociationPoints = _r.Count;
            }
            //we assume the input array is well prepared and we don't check for errors
            if (this._NumOfDissociationPoints != -1 && _r.Count != this._NumOfDissociationPoints)
            {
                throw new Exception("array are not set correctly with not correct number of elements");
            }
            _RUDissociation = _r.ToArray();//it is not very important to be in order ???
        }

        public void SetConc(double _c)
        {
            this._Conc=_c;
        }

        public double[] AssociationTimes
        {
            get { return this._TimeAssociation; }
        }
        public double[] DissociationTimes
        {
            get { return this._TimeDissociation; }
        }
        public double[] AssociationRUs
        {
            get { return this._RUAssociation; }
            //set { this._RUAssociation = value; }
        }
        public double[] DissociationRUs
        {
            get { return this._RUDissociation; }
            //set { this._RUDissociation = value; }
        }
        public double Concentration
        {
            get { return this._Conc; }
        }
        private double _Conc;
        private int _NumOfAssociationPoints;
        private int _NumOfDissociationPoints;
        private double [] _TimeAssociation;
        private double [] _TimeDissociation;
        private double [] _RUAssociation;
        private double [] _RUDissociation;
    }
    
}
