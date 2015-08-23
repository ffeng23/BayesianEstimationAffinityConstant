using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Models;
namespace TestingGibbsSamplerOnQuadraticEquation
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Hello world!!!");

            //Start doing the 
            FitController fc = new QuadraticFitController();

            fc.SetupModel();
            List<List<double>> dist=fc.Run(1);

            Console.WriteLine("Start doing the fitting of ELISTAT");
            FitController fc_ELISTAT = new ELISTATFitController();

            fc_ELISTAT.SetupModel();
            List<List<double>> dist_ELISTAT = fc_ELISTAT.Run(1);

            Console.WriteLine("Start doing the fitting of ELISTAT Full Model");
            FitController fc_ELISTAT_QM = new ELISTATQuadraticFitController();

            fc_ELISTAT_QM.SetupModel();
            List<List<double>> dist_ELISTAT_QM = fc_ELISTAT_QM.Run(500);

            Console.WriteLine("\n\nPress ENTER to quit.............");
            Console.ReadLine();
        }
    }
}
