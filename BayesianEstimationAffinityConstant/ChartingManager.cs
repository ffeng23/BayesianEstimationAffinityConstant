using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows.Forms.DataVisualization.Charting;

using Meta.Numerics.Matrices;
using System.Windows.Forms;
using System.Drawing;


namespace BayesianEstimationAffinityConstant
{
    public class ChartingManager
    {
        public ChartingManager()
        {
            //build color table for later use
            colorTable = new Dictionary<int, Color>();
            colorTable.Add(0, Color.Red);
            colorTable.Add(1, Color.Green);
            colorTable.Add(2, Color.Blue);
            colorTable.Add(3, Color.Brown);
            colorTable.Add(4, Color.SeaGreen);
        }
        //***********charting for cell division
        public void DrawTracePlots(List<List<double>> _xData, List<List<double>> _yData, List<string> _title, List<string> _xlab, List<string> _ylab, bool drawLine = false)
        {
            foreach (Control c in pChart.Controls)
                pChart.Controls.Remove(c);
            cChart = new Chart();
            for (int i = 0; i < _xData.Count; i++)
            {

                ChartArea chartArear1 = new ChartArea(""+i);
                cChart.ChartAreas.Add(chartArear1);
                //double[] x; double[] y;
                //x = time.ToArray();
                //y = cellNumber.ToArray();
                drawTracePlot(_xData[i], _yData[i], chartArear1, _title[i], _xlab[i], _ylab[i], "", 0, drawLine);
                chartArear1.Position.X =1;
                chartArear1.Position.Y = 5 + i *16;
                chartArear1.Position.Width = 99;
                chartArear1.Position.Height = 16F;
                if (i > 0)
                {
                    chartArear1.AlignWithChartArea = "0";
                    chartArear1.AlignmentOrientation = AreaAlignmentOrientations.Vertical;
                    chartArear1.AlignmentStyle = AreaAlignmentStyles.All;
                }

               
            }

            //*********************

            cChart.Titles.Add("Title_1");
            cChart.BackColor = Color.White;
            

            cChart.Titles[0].Text = "Trace Plots";
            // Set chart control location
            //cChart.Location = new System.Drawing.Point(1,1);

            // Set Chart control size
            cChart.Size = new System.Drawing.Size( 550,810 );
            pChart.Controls.Add(cChart);
            
            return;
        }

        private void drawTracePlot(List<double> x, List<double> y, ChartArea cA, string title, string xlab, string ylab,
            string seriesName = "", int _color = 0, bool drawLine = false)
        {
            if (x.Count() == 0 && y.Count() == 0)
            {
                //no point to draw, return
                return;
            }
            //draw points first
            int n = x.Count() <= y.Count() ? x.Count() : y.Count();
            Series s = new Series(seriesName);

            cChart.Series.Add(s);
            for (int i = 0; i < n; i++)
            {
                s.Points.AddXY(x[i], y[i]);
            }

            // Add series to the chart
            if (drawLine)
            {
                s.ChartType = SeriesChartType.FastLine;
            }
            else
            {
                s.ChartType = SeriesChartType.Point;
            }
            s.ChartArea = cA.Name;
            s.MarkerSize = 4;
            s.MarkerStyle = MarkerStyle.Circle;
            s.Color = colorTable[_color % colorTable.Count];
            //string[] chartAxisLabel = { "x", "y", "z" };
            cA.AxisY.Minimum = y.Min() ;
            cA.AxisY.Maximum = y.Max();
            cA.AxisX.Minimum = x.Min() ;
            cA.AxisX.Maximum = x.Max() ;
            cA.AxisY.Interval = (cA.AxisY.Maximum - cA.AxisY.Minimum) / 3;
            cA.AxisX.Title = xlab;
            cA.AxisY.Title = ylab;
            cA.AxisX.TitleFont = new Font("Arial", 8);
            cA.AxisY.TitleFont = new Font("Arial", 8);

            cA.AxisX.LabelStyle.Font = new Font("Arial", 7);
            cA.AxisY.LabelStyle.Font = new Font("Arial", 7);

            cA.AxisX.LabelStyle.Format = "##,#";
            cA.AxisY.LabelStyle.Format = "G2";
            cA.Visible = true;
            
            //cA.Position.Auto = true;
            //cA.Position.X = 3 ;
            //cA.Position.Y = 10 ;
            //cA.Position.Height = 90;//80-offset*5;
            //cA.Position.Width = 90;// 80 - offset * 5;
            //Color[] cls = { Color.LightCyan, Color.LightCyan, Color.LightCyan };
            //cA.BackColor = Color.Red;
            cA.AxisX.MajorGrid.Enabled=false;
            cA.AxisY.MajorGrid.Enabled = false ;
            cA.BorderWidth =1;
            cA.BorderDashStyle = ChartDashStyle.Solid ;
            cA.BorderColor = Color.Black;
        }
        


        /// <summary>
        /// setting the charting control hold on the wpf
        /// </summary>
        /// <param name="p">the panel hold where to draw the chart</param>
        public void setChartingPanel(Panel p)
        {
            pChart = p;

        }

        /// <summary>
        /// access drawing panel for charting error
        /// </summary>
        public Panel PChart
        {
            get { return pChart; }

        }

        private Chart cChart;
        private Panel pChart;
        private Dictionary<int, Color> colorTable;
    }//end of class
}
