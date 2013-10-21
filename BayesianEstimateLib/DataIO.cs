using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Text.RegularExpressions;

namespace BayesianEstimateLib
{
    /// <summary>
    /// this is the class taking care of data input/output
    /// </summary>
    public class DataIO
    {
        private const int BUFFER_LENGTH = 5000;
        public DataIO()
        {
        }
        /// <summary>
        /// 
        /// this assumes a fixed formation, it arranges the data in columns, first column is time data, then follows a variable number
        /// of data columns (RU columns). It has a head, but is skipped. the output is the a dictionary, keyed by integer, first one, key=0 
        /// is the time, then the RUs, keyed by 1,2,.....
        /// 
        /// </summary>
        /// <param name="fileName">string filename cotaining the diectory path</param>
        /// <returns></returns>
        public static Dictionary<int, List<Double>> ReadDataTable(String fileName)
        {
            string line;
            int nRow = 0;
            Dictionary<int, List<Double>> output = new Dictionary<int, List<double>>();
            FileInfo src = new FileInfo(fileName);
            if (!src.Exists)
            {
                throw new System.Exception("file doesn't exist");
               // Environment.Exit(-1);
                
            }
            TextReader rin = src.OpenText();
            line = rin.ReadLine();
            string[] buf = new string[BUFFER_LENGTH];
            string tems;
            //need to use the first header row to determine the number of columns
            int nCol = Regex.Matches(line, "\t").Count + 1;
            //here we assume at least the header column is a good format in terms of number of columns
            for (int i = 0; i < nCol; i++)
            {
                output.Add(i, new List<Double>());
            }

            line = rin.ReadLine();
            while (line != null)
            {
                nRow++;
                int runningColNum = Regex.Matches(line, "\t").Count + 1;
                if (runningColNum != nCol)
                {
                    throw new System.Exception("wrong file format, the colnum numbers are variable!");
                    //Environment.Exit(-1);
                }
                buf = line.Split('\t');
                List<string> temp = new List<string>();
                for (int i = 0; i < nCol; i++)
                {
                    tems = Regex.Replace(buf[i], "^\\s+", "");
                    tems = Regex.Replace(tems, "\\s+$", "");
                    if (tems.Length != 0)
                        output[i].Add(Convert.ToDouble(tems));
                    else //add nan
                        output[i].Add(Double.NaN);
                }
                //element.Add(temp);
                line = rin.ReadLine();
            }
            rin.Close();
            return output;
        }
        /// <summary>
        /// read the datatable and return a dictionary keyed by the header (if no header, then keyed by interger indiciation the index of the column
        /// </summary>
        /// <param name="_fileName"></param>
        /// <param name="_Header"></param>
        /// <param name="_SepChar"></param>
        /// <param name="_SkipLine"></param>
        /// <returns></returns>
        public static Dictionary<string, List<Double>> ReadDataTable(String _fileName, bool _Header, Char _SepChar, int _SkipLine=0)
        {
            string line;
            int nRow = 0;
            Dictionary<string, List<Double>> output = new Dictionary<string, List<double>>();
            FileInfo src = new FileInfo(_fileName);
            if (!src.Exists)
            {
                throw new System.Exception("file doesn't exist");
                // Environment.Exit(-1);

            }
            TextReader rin = src.OpenText();
            
            //first skip the _skiplines

            do
            {
                if (nRow >= _SkipLine)
                    break;
                line = rin.ReadLine();
                nRow++;
            } while (line != null);

            //need to use the first header row to determine the number of columns
            line = rin.ReadLine();
            string[] buf = new string[BUFFER_LENGTH];
            string tems;
            int nCol = Regex.Matches(line, _SepChar.ToString()).Count + 1;
            if (line == null)
            {
                return null;
            }
            buf = line.Split(_SepChar);
            if (!_Header) //no header
            {
                //here we assume at least the header column is a good format in terms of number of columns
                for (int i = 0; i < nCol; i++)
                {
                    output.Add(i.ToString(), new List<Double>());
                }
            }
            else //read in the header
            {
                for (int i = 0; i < nCol; i++)
                {
                    output.Add(buf[i], new List<Double>());
                }
            }
            line = rin.ReadLine();
            while (line != null)
            {
                nRow++;
                int runningColNum = Regex.Matches(line, _SepChar.ToString()).Count + 1;
                if (runningColNum != nCol)
                {
                    throw new System.Exception("wrong file format, the colnum numbers are variable!");
                    //Environment.Exit(-1);
                }
                buf = line.Split(_SepChar);
                List<string> temp = new List<string>();
                for (int i = 0; i < nCol; i++)
                {
                    tems = Regex.Replace(buf[i], "^\\s+", "");
                    tems = Regex.Replace(tems, "\\s+$", "");
                    if (tems.Length != 0)
                        output[output.Keys.ToList()[i]].Add(Convert.ToDouble(tems));
                    else //add nan
                        output[output.Keys.ToList()[i]].Add(Double.NaN);
                }
                //element.Add(temp);
                line = rin.ReadLine();
            }
            rin.Close();
            return output;
        }

        public static void WriteDataTable(String _filename, Dictionary<int, List<double>> _lst)
        {
            StreamWriter writer = new StreamWriter(_filename);
            List<int> header = _lst.Keys.ToList();
            for (int i = 0; i < header.Count; i++)
            {
                writer.Write(header[i]);
                   if(i!=header.Count-1)
                       writer.Write("\t");
            }
            writer.Write("\r\n");
            for (int i = 0; i < _lst[0].Count; i++)
            {
                for (int j = 0; j < header.Count; j++)
                {
                    writer.Write(_lst[header[j]][i]);
                    if (j != header.Count - 1)
                        writer.Write("\t");
                }
                writer.Write("\r\n");
            }
            
            writer.Close();
        }
        public static void WriteDataTable(List<double> lst1, List<double> lst2, string _filename, List<string> header)
        {
            StreamWriter writer = new StreamWriter(_filename);
            writer.WriteLine(header[0] + "\t" + header[1]);
            for (int i = 0; i < lst1.Count; i++)
            {
                writer.WriteLine(lst1[i] + "\t" + lst2[i]);
            }
            writer.Close();

        }

        public static void Impute(List<Double> _array)
        {
            for (int i = 0; i < _array.Count; i++)
            {
                if(double.IsNaN(_array[i]))
                {
                    if(i==0)
                        _array[i]=0;
                    else //need to compute the average between the one around it;
                    {
                        if(i==_array.Count-1) //last one
                        {
                            if(_array.Count==2)
                            {
                                _array[i]=_array[i-1];
                            }
                            else
                            {
                                _array[i]=_array[i-1]+(_array[i-1]-_array[i-2]);
                            }
                        }
                        else //this is not the last one
                        {
                            double lower=_array[i-1];
                            int runningIndex=1;
                            while(double.IsNaN(_array[i+runningIndex])&&i+runningIndex<_array.Count)
                            {
                                runningIndex++;
                            }

                            if(i+runningIndex>=_array.Count)//this is a sick situation, we are having all trailing NaNs.
                            {
                                _array[i]=_array[i-1]+(_array[i-1]-_array[i-2]);
                            }
                            else //we are finding some good double nubmer follwing this current one
                            {
                                _array[i]=1/(runningIndex+1)*(_array[i+runningIndex]-_array[i-1])+_array[i-1];
                            }

                         }//end of else

                    }
                }
            }
        }
    }//end of the class
}
