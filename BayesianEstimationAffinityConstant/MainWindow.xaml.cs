using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

using System.IO;
using BayesianEstimateLib;
using System.Windows.Forms;
using Models;

namespace BayesianEstimationAffinityConstant
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();

            chtManager = new ChartingManager();
        }
        public System.Windows.Controls.TextBox LogTextBlock
        {
            get { return this.tBk_Log; }
        }

        private NumericalIntegrationOfDynamics nid;
        private double ka;
        private double kd;
        private double conc;
        private double Rmax;
        private double R0;
        private double durationA;
        private double durationD;
        private double kM;
        private double deltaT;
        //private double Rmax;

        private List<double> _time;
        private List<double> _ru;

        private double var;

        private void btnRunSim_Click(object sender, RoutedEventArgs e)
        {
            //getting values first
            ka = Convert.ToDouble(Tbx_ka.Text);
            kd = Convert.ToDouble(Tbx_kd.Text);
            conc = Convert.ToDouble(Tbx_Conc.Text);
            Rmax = Convert.ToDouble(Tbx_Rmax.Text);
            R0 = Convert.ToDouble(Tbx_R0.Text);
            //durationA = 400;// Convert.ToDouble(Tbx_DurationA.Text);
            //durationD = 400;// Convert.ToDouble(Tbx_DurationD.Text);
            kM  = Convert.ToDouble(Tbx_kM.Text);
            deltaT = Convert.ToDouble(Tbx_deltaT.Text);
            var = Convert.ToDouble(Tbx_Var.Text); 
            //run simulation
            nid = new NumericalIntegrationOfDynamics(ka, kd, conc, Rmax, R0, durationA, durationD, kM, deltaT);
            double r0 = 0;
            List<string> header = new List<string>(2);
            header.Add("time");
            header.Add("RU");
            //doing the simulation for different Rmax

            nid.run_Attach();
            nid.R0 = nid.RU_Attach[nid.RU_Attach.Count - 1];
            nid.run_Detach();


            Console.WriteLine("write output..................");
            LogTextBlock.Text += "\nDoing simulation........";
            writeOutput(nid.Time_Attach, nid.RU_Attach, "simulation_attach.txt", header);
            writeOutput(nid.Time_Detach, nid.RU_Detach, "simulation_detach.txt", header);


            //writeOutput(nid.Time_Attach, nid.RU_Attach, "simulation_attach_Rmax"+(_Rmax-stepOfRmax*i) +".txt", header);
            //writeOutput(nid.Time_Detach, nid.RU_Detach, "simulation_detach_Rmax" + (_Rmax - stepOfRmax*i) + ".txt", header);

            Console.WriteLine("add noise...............");
            LogTextBlock.Text +="add noise.............";
            r0 = nid.RU_Attach[nid.RU_Attach.Count - 1];
            nid.addNoise(Math.Sqrt(var));

            writeOutput(nid.Time_Attach, nid.RU_Attach, "simulation_attach_noise.txt", header);
            writeOutput(nid.Time_Detach, nid.RU_Detach, "simulation_detach_noise.txt", header);

            Console.WriteLine("done............");

            //now we need to shorten the data input
            List<double> shortT_A = new List<double>();
            List<double> shortR_A = new List<double>();
            List<double> shortT_D = new List<double>();
            List<double> shortR_D = new List<double>();
            for (int i = 0; i < nid.Time_Attach.Count; i = i + 3)
            {
                shortT_A.Add(nid.Time_Attach[i]);
                shortR_A.Add(nid.RU_Attach[i]);
            }
            for (int i = 0; i < nid.Time_Detach.Count; i = i + 3)
            {
                shortT_D.Add(nid.Time_Detach[i]);
                shortR_D.Add(nid.RU_Detach[i]);
            }
            
            writeOutput(shortT_A, shortR_A, "simulation_attach_noiseShort.txt", header);
            writeOutput(shortT_D, shortR_D, "simulation_detach_noiseShort.txt", header);

            LogTextBlock.Text += "\rdone.....";
            Console.WriteLine("done..........");
        }

        private void writeOutput(List<double> lst1, List<double> lst2, string _filename, List<string> header)
        {
            StreamWriter writer = new StreamWriter(_filename);
            writer.WriteLine(header[0]+"\t"+header[1]);
            for (int i = 0; i < lst1.Count; i++)
            {
                writer.WriteLine(lst1[i]  +"\t" + lst2[i]);
            }
            writer.Close();

        }

        private void BtnRunBayesian_Click(object sender, RoutedEventArgs e)
        {
            this.WindowState = WindowState.Minimized;
            this.Topmost = true;
            this.Topmost = false;
            //this.IsVisible = false;
            Console.WriteLine("*****testing. fitting the model......");

            //look to see whether we want to use the new simulated input or the one from the file
            string fileAttach = "simulation_attach_noiseShort.txt";
            string fileDetach = "simulation_detach_noiseShort.txt";
            
            Console.WriteLine("*****reading the input from disk*****");
            
            fileAttach = Tbx_AttachDataFile.Text ;
            fileDetach = Tbx_DetachDataFile.Text;
            //we need to read in the input
            Dictionary<string, List<double>> inputDataShort = BayesianEstimateLib.DataIO.ReadDataTable(fileAttach, true, '\t', 0);
            List<double> shortT_A = inputDataShort["time"];
            List<double> shortR_A = inputDataShort["RU"];
            
            //Console.WriteLine("reading \"" + args[1] + "\"........");
            inputDataShort = BayesianEstimateLib.DataIO.ReadDataTable(fileDetach, true, '\t', 0);
            List<double> shortT_D = inputDataShort["time"];
            List<double> shortR_D = inputDataShort["RU"];
            Console.WriteLine("Done");
            //}
            SprModel sprm = new SprModel(shortT_A, shortR_A, shortT_D, shortR_D);
            //getting values first
            ka = Convert.ToDouble(Tbx_ka.Text);
            kd = Convert.ToDouble(Tbx_kd.Text);
            conc = Convert.ToDouble(Tbx_Conc.Text);
            Rmax = Convert.ToDouble(Tbx_Rmax.Text);
            R0 = Convert.ToDouble(Tbx_R0.Text);
            //durationA = Convert.ToDouble(Tbx_DurationA.Text);
            //durationD = Convert.ToDouble(Tbx_DurationD.Text);
            kM = Convert.ToDouble(Tbx_kM.Text);
            deltaT = Convert.ToDouble(Tbx_deltaT.Text);
            var = Convert.ToDouble(Tbx_Var.Text); 

            sprm.SetAllParamters(ka, kd, kM, conc, Rmax, R0, var);

            List<int> lstFunc = new List<int>();
            lstFunc.Add(0/*ka*/);
            lstFunc.Add(1/*kd*/);
            lstFunc.Add(3/*conc*/);
            lstFunc.Add(4/*Rmax*/);
            lstFunc.Add(5);/*R0*/
            lstFunc.Add(6/*var*/);
            sprm.setFunctionDelegateForUpdating(lstFunc);

            Console.WriteLine("setting up the lower/upper bounds of parameters..........");
            List<List<double>> bounds = new List<List<double>>();
            bounds.Add(new List<double> { 1E5, 1E9 });//ka, (1, 1E16)
            bounds.Add(new List<double> { 0, 2E0 });//kd, (0, 100)
            //bounds.Add(new List<double> { 1E2, 1E10 });//kM, (0, 1E15)
            bounds.Add(new List<double> { 0, 1E-4 });//conc, (0, 10)
            bounds.Add(new List<double> { 5E0, 5E2 });//Rmax, (0, 1E5)
            bounds.Add(new List<double> { 1, 1E3 });//R0, (0, 1E5)
            bounds.Add(new List<double> { 0, 20 });//var, (0, 3000)

            Console.WriteLine("building the gibbs sampler...........");
            GibbsSampler.GibbsSampler gb = new GibbsSampler.GibbsSampler(
                /*init values*/new List<double> { Convert.ToDouble(Tbx_ka.Text), Convert.ToDouble(Tbx_kd.Text), Convert.ToDouble(Tbx_Conc.Text),
                   Convert.ToDouble(Tbx_Rmax.Text), Convert.ToDouble(Tbx_R0.Text), Convert.ToDouble(Tbx_Var.Text)}
                   , sprm.updateFunctionDistribution, bounds);
            //GibbsSampler.GibbsSampler gb = new GibbsSampler.GibbsSampler(new List<double> {1545903.57844479,	0.000322900628507969,	9.91536866021184E-09,	29.8920610679908,	29.5843847101981,	1.65389525923268}, sprm.updateFunctionDistribution, bounds);
            Console.WriteLine("now ready to run estimation........");
            //List<List<double>> output = new List<List<double>>();
            runResultData  = gb.Run(Convert.ToInt32(Tbx_TotalSteps.Text));
            

            Console.WriteLine("writing the output file.......");
            StreamWriter writer = new StreamWriter("learReg.txt");
            //writer.WriteLine("ka\tkb\tkM\tconc\tRmax\tR0\tVar");
            for (int i = 0; i < runResultData[0].Count; i++)
            {
                for (int j = 0; j < runResultData.Count; j++)
                {
                    writer.Write(runResultData[j][i]);
                    if (j != runResultData.Count - 1)
                    {
                        writer.Write("\t");
                    }
                    else
                        writer.Write("\r\n");
                }
            }
            Console.WriteLine("done!!!!!!!");
            writer.Close();

            RenderChartingTracePlots(runResultData , pChartTracePlot);

            Console.WriteLine("Done...........");

            this.WindowState = WindowState.Normal;
            this.Topmost = false;
            this.Topmost = true;
            this.Focus();

            


        }

        private void Btn_AttachingFile_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog();

            //int size = -1;
            
            openFileDialog1.DefaultExt = ".txt";
            openFileDialog1.Filter = "Text documents (.txt)|*.txt";
 
            // Display OpenFileDialog by calling ShowDialog method
            Nullable<bool> result = openFileDialog1.ShowDialog();
 
            // Get the selected file name and display in a TextBox
            if (result == true)
            {
                // Open document
                string filename = openFileDialog1.FileName;
                Tbx_AttachDataFile.Text = filename;
             }
            //Console.WriteLine(size); // <-- Shows file size in debugging mode.
            //Console.WriteLine(result); // <-
        }

        private void Btn_DetachingFile_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog();

            //int size = -1;

            openFileDialog1.DefaultExt = ".txt";
            openFileDialog1.Filter = "Text documents (.txt)|*.txt";

            // Display OpenFileDialog by calling ShowDialog method
            Nullable<bool> result = openFileDialog1.ShowDialog();

            // Get the selected file name and display in a TextBox
            if (result == true)
            {
                // Open document
                string filename = openFileDialog1.FileName;
                Tbx_DetachDataFile.Text = filename;
            }

        }

        private void Expdr_Expanded(object sender, RoutedEventArgs e)
        {
            string text="";
            if (sender == Expdr_ka)
            {
                text = "ka -- (1/(Ms))   ";
            }
            if (sender == Expdr_kd)
            {
                text = "kd -- (1/s)   ";
            }
            if (sender == Expdr_A0)
            {
                text = "A0 -- (M)   ";
            }
            if (sender == Expdr_Rmax)
            {
                text = "Rmax -- (RU)   ";
            }
            if (sender == Expdr_R0)
            {
                text = "R0 -- (RU)   ";
            }
            if (sender == Expdr_kM)
            {
                text = "kM -- (RU/(Ms))   ";
            }
            ((Expander)sender).Header =  text;
            ((Expander)sender).Background = Brushes.LightGray;

            System.Windows.Controls.Panel.SetZIndex(((UIElement)((Expander)sender).Content), 1000); Topmost = true;
        }

        private void Expdr_Collapsed(object sender, RoutedEventArgs e)
        {
            ((Expander)sender).Header = "";
            ((Expander)sender).Background = null;
        }

        private void Tbx_kM_TextChanged(object sender, TextChangedEventArgs e)
        {
            UpdateDiffusionAndkM();   
        }
        private void UpdateDiffusionAndkM()
        {
            if (Ckb_kM != null && Ckb_kM.IsChecked == true)
            {
                //in this case, we don't update anything
                return;
            }
            //Console.WriteLine("in here");
            //take input
            double temp = 0;
            double mw = 0;
            double L1 = 0;
            double L2 = 0;
            double h = 0;
            double w = 0;
            double FlowRate = 0;
            double D = 0;
            bool nullFlag = false;
            try
            {
                if (this.Tbx_kM_Temp == null)
                {
                    nullFlag = true;
                }
                else
                {
                    temp = Convert.ToDouble(this.Tbx_kM_Temp.Text);
                }
                if (this.Tbx_kM_MW == null)
                {
                    nullFlag = true;
                }
                else
                {
                    mw = Convert.ToDouble(this.Tbx_kM_MW.Text);
                }

                if (this.Tbx_kM_l1 == null)
                {
                    nullFlag = true;
                }
                else
                {
                    L1 = Convert.ToDouble(this.Tbx_kM_l1.Text);
                }

                if (this.Tbx_kM_l2 == null)
                {
                    nullFlag = true;
                }
                else
                {
                    L2 = Convert.ToDouble(this.Tbx_kM_l2.Text);
                }
                if (this.Tbx_kM_h == null)
                {
                    nullFlag = true;
                }
                else
                {
                    h = Convert.ToDouble(this.Tbx_kM_h.Text);
                }

                if (this.Tbx_kM_w == null)
                {
                    nullFlag = true;
                }
                else
                {
                    w = Convert.ToDouble(this.Tbx_kM_w.Text);
                }
                if (this.Tbx_kM_flowRate == null)
                {
                    nullFlag = true;
                }
                else
                {
                    FlowRate = Convert.ToDouble(this.Tbx_kM_flowRate.Text);
                }
            }
            catch (Exception ex)
            {
                //throw (ex);
                //Console.WriteLine(ex.ToString());
                Console.WriteLine("*******Error: invalid input format!!!!!!");
                Console.Write("\a");
                nullFlag = true;
            }
            if (nullFlag)
                return;

            // System.Windows.MessageBox.Show("Yes");

            if (this.CkB_DiffusionCoefficient.IsChecked == true)
            {//we need to take value as user input
                D = Convert.ToDouble(this.Tbx_kM_D.Text);
            }
            else
            {
                D = MassTransportCoefficient_kM.CalculateDiffusionCoefficient(temp, mw);
                Tbx_kM_D.Text = "" + D.ToString("e2");
            }
            double kM = MassTransportCoefficient_kM.CalculatekM(temp, mw, FlowRate, L1, L2, h, w, D);
            double kM_prime = MassTransportCoefficient_kM.CalculatekM_prime(mw, kM);
            this.Tbx_kM.Text = "" + kM_prime.ToString("e2");
        }

        private void CkB_DiffusionCoefficient_Checked(object sender, RoutedEventArgs e)
        {
            if (this.CkB_DiffusionCoefficient.IsChecked == true)
            {
                //we will allow the take diffusion coefficient as input
                this.Tbx_kM_D.IsEnabled = true;

            }
            else
            {
                this.Tbx_kM_D.IsEnabled = false;
            }
        }

        private void Ckb_kM_Checked(object sender, RoutedEventArgs e)
        {
            Tbx_kM.IsEnabled = true;
        }

        private void Ckb_kM_Unchecked(object sender, RoutedEventArgs e)
        {
            Tbx_kM.IsEnabled = false;

            UpdateDiffusionAndkM();
        }
        //**************starting doing the charting for traceplot
        private void RenderChartingTracePlots(List<List<double>> _outdata, System.Windows.Forms.Panel p)
        {
            if (chtManager == null)
                return;
            chtManager.setChartingPanel(p);
            int burnInSteps=Convert.ToInt32(this.Tbx_BurnIn.Text);
            //need to prepare the output data for plotting
            List<List<double>> xdata = new List<List<double>>();
            List<string> title = new List<string>(_outdata.Count) { "", "", "", "", "", "" }; //{ "ka", "kd", "conc", "Rmax", "R0", "Var" };
            List<string> xlab=new List<string>{"", "", "", "", "", ""};//{"ka", "kd", "conc", "Rmax", "R0", "Var"};
            List<string> ylab=new List<string>{"ka", "kd", "conc", "Rmax", "R0", "Var"};//{"", "", "", "", "", ""};
            List<List<double>> burnInData = new List<List<double>>(xdata.Count);

            for (int i = 0; i < _outdata.Count; i++)
            {
                List<double> xdata_x = new List<double>(_outdata[i].Count);
                List<double> burnIndata_y = new List<double>();
                for (int j = 0; j < _outdata[i].Count; j++)
                {
                    if (j >= burnInSteps)
                    {
                        xdata_x.Add(j + 1);
                        burnIndata_y.Add(_outdata[i][j]);
                    }
                }
                xdata.Add(xdata_x);
                burnInData.Add(burnIndata_y);
                
            }

            if (burnInData[0].Count <= 1)
            {
                return;
            }
            chtManager.DrawTracePlots(xdata,_outdata, title, xlab, ylab,true );
            printDistribution(burnInData);

        }
        private void printDistribution(List<List<double>> _data)
        {
            //need to write the summary table on the blockd
            List<double> b_mean=new List<double>(_data.Count );
            List<double> b_lower5Percentile = new List<double>(_data.Count);
            List<double> b_upper5Percentile = new List<double>(_data.Count);
            this.tBk_Log.Text += "\n*******************************************\n";
            this.tBk_Log.Text += "Bayesian Estimation Summary\n";
            this.tBk_Log.Text += "=======================================\n";
            this.tBk_Log.Text += " \tmean\t5% Percentile\t95% Percentile\n";
            this.tBk_Log.Text += "------------------------------------------------------\n";
            

            //now get them done

            for (int i = 0; i < _data.Count; i++)
            {
                switch (i)

                {
                    case 0:

                        this.tBk_Log.Text += "ka\t";
                        break;
                    case 1:
                        this.tBk_Log.Text += "kd\t";
                        break;
                    case 2:
                        this.tBk_Log.Text += "Conc\t";
                        break;
                    case 3:
                        this.tBk_Log.Text += "Rmax\t";
                        break;
                    case 4:
                        this.tBk_Log.Text += "R0\t";
                        break;
                    case 5:
                        this.tBk_Log.Text += "Var\t";
                        break;
                    default:
                        Console.Beep();
                        Console.Write("unidentified paramter value!!!!!!");
                        
                        //Environment.exit(-1);
                        break;
                }
                b_mean.Add( 0);
                for (int j = 0; j < _data[i].Count; j++)
                {
                    b_mean[i] += _data[i][j];
                }
                b_mean[i] /= _data[i].Count;
                _data[i].Sort();
                int index=((int) (_data[i].Count * 0.05));
                b_lower5Percentile.Add(_data[i][index]);
                index = ((int)(_data[i].Count * 0.95));
                b_upper5Percentile.Add(_data[i][index]);
                this.tBk_Log.Text += b_mean[i].ToString("e4") + "\t" + b_lower5Percentile[i].ToString("e4") + "\t" + b_upper5Percentile[i].ToString("e4") + "\n";
                
                
            }
            //write to the log area
            this.tBk_Log.Text += "===================================\n";
            this.tBk_Log.Text += "Total steps: " + this.runResultData[0].Count + ";\tBurnIn steps: " + (this.runResultData[0].Count- _data[0].Count)+".\n\n";

        }

        private ChartingManager chtManager;
        private List<List<double>> runResultData;

        private void Tbx_BurnIn_TextChanged(object sender, TextChangedEventArgs e)
        {
            if (pChartTracePlot == null)
                return;
            //redraw the trace plot
            if (runResultData ==null ||runResultData.Count ==0)
            {
                return;
            }
            RenderChartingTracePlots(runResultData, pChartTracePlot);
        }
    }//end of class
}//end of namespace
