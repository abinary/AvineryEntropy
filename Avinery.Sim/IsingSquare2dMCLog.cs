using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Avinery.Sim
{
    public class IsingSquare2dMCLog
    {
        public IsingSquare2dMCLog(int numOfEntries, int dn_N = 4)
        {
            NumOfEntries = numOfEntries;
            Mean = new double[numOfEntries];
            Variance = new double[numOfEntries];
            Trajectories = new int[numOfEntries][,];

            B1 = new int[1 << 2];
            B2 = new int[1 << 4];
            B3 = new int[1 << 6];
            B4 = new int[1 << 8];
            B5 = new int[1 << 10];
            B6 = new int[1 << 12];
            B7 = new int[1 << 14];

            B1x = new int[1 << 2];
            B2x = new int[1 << 4];
            B3x = new int[1 << 6];
            B4x = new int[1 << 8];
            B5x = new int[1 << 10];
            B6x = new int[1 << 12];
            B7x = new int[1 << 14];

            this.dn_N = dn_N;
            //dn_N = 3;
            dn_H = new int[1 << 2* dn_N];
            dn_H_excluding_O = new int[1 << 2 * dn_N - 1];
        }

        public int NumOfEntries { get; set; }

        public double[] Mean { get; set; }
        public double[] Variance { get; set; }

        public readonly int dn_N = 4;
        public int[] dn_H { get; set; }
        public int[] dn_H_excluding_O { get; set; }
        public double dn_Entropy;

        public int[] B1 { get; set; }
        public int[] B2 { get; set; }
        public int[] B3 { get; set; }
        public int[] B4 { get; set; }
        public int[] B5 { get; set; }
        public int[] B6 { get; set; }
        public int[] B7 { get; set; }

        public int[] B1x { get; set; }
        public int[] B2x { get; set; }
        public int[] B3x { get; set; }
        public int[] B4x { get; set; }
        public int[] B5x { get; set; }
        public int[] B6x { get; set; }
        public int[] B7x { get; set; }

        public int[][,] Trajectories { get; set; }
        public int[,] SampledSpins { get; set; }

        public void FinalizeCalculations()
        {
            double dn_H_sum = dn_H.Sum();
            double dn_H_excluding_O_sum = dn_H_excluding_O.Sum();

            // These are Log2 operations
            var s2 = dn_H.Select(x => x/dn_H_sum).Select(x => x*Math.Log(x + 1e-15, 2.0)).Sum();
            var s1 = dn_H_excluding_O.Select(x => x / dn_H_excluding_O_sum).Select(x => x * Math.Log(x + 1e-15, 2.0)).Sum();

            this.dn_Entropy = s2 - s1;
        }
    }
}
