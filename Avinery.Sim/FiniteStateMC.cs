using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;

namespace Avinery.Sim
{
    public class FiniteStateMCLog
    {
        public FiniteStateMCLog(int numOfEntries)
        {
            NumOfEntries = numOfEntries;
            VisitLog = new int[numOfEntries];
        }

        public int NumOfEntries { get; set; }
        public int[] VisitLog { get; set; }

        public double[] MeanOccupancy { get; set; }
        public double Entropy { get; set; }
        public double TheoreticalEntropy { get; set; }
    }

    public class FiniteStateMC : MCBase
    {
        public FiniteStateMC()
        {
            T = 1.0;
        }

        private int[] _stateVisits;
        private double[] _stateEnergies;

        public int CurrentState { get; set; }

        public double T { get; set; }

        public double[] StateEnergies
        {
            get { return _stateEnergies; }
            set { _stateEnergies = value; ResetStateVisits(); }
        }

        public int[] StateVisits
        {
            get { return (int[])_stateVisits.Clone(); }
        }

        public void ResetStateVisits()
        {
             _stateVisits = new int[_stateEnergies.Length];
        }

        public int OutputStyle { get; set; }
        public string OutputFileName { get; set; }
        public Stream OutputStream { get; set; }

        public FiniteStateMCLog RandomVisits(int numOfSteps)
        {
            var log = new FiniteStateMCLog(numOfSteps);

            var inverseTemp = 1/T;
            var stateProbabilities = _stateEnergies.Select(e => Math.Exp(-(e* inverseTemp))).ToArray();

            var p = StateEnergies.Divide(T).Negative().Exp();
            p = p.Divide(p.Sum());

            var cumulativeProbability = p.CumSum();

            Stream trajectoriesStream = null;

            if (OutputStream != null)
                trajectoriesStream = OutputStream;
            else if (OutputFileName != null)
                trajectoriesStream = new FileStream(OutputFileName, FileMode.Create, FileAccess.Write, FileShare.Read);

            Action<int> writeFunction;
            BitWriter bs = null;
            var numOfBits = (int)Math.Ceiling(Math.Log(StateEnergies.Length, 2.0));

            switch (OutputStyle)
            {
                case 1:
                    writeFunction = (state => { trajectoriesStream.WriteByte((byte) (1 << state)); });
                    break;

                case 2:
                    bs = new BitWriter(trajectoriesStream);
                    writeFunction = (state => { bs.WriteBits((ulong)state, numOfBits); });
                    break;

                default:
                    writeFunction = (state => { trajectoriesStream.WriteByte((byte)state); });
                    break;
            }

            for (int stepIndex = 0; stepIndex < numOfSteps; stepIndex++)
            {
                var randomNum = _random.NextDouble();

                CurrentState = cumulativeProbability.IndexOfFirstGreater(randomNum);

                _stateVisits[CurrentState]++;
                log.VisitLog[stepIndex] = CurrentState;

                // Write trajectories
                if (trajectoriesStream != null)
                    writeFunction(CurrentState);
            }

            // Dispose of the stream only if this funcion created it
            if (trajectoriesStream != null && OutputStream == null)
                trajectoriesStream.Dispose();

            double sum = _stateVisits.Sum();
            log.MeanOccupancy = _stateVisits.Select(v => v / sum).ToArray();
            log.Entropy = -log.MeanOccupancy.Select(m => (m + 1e-16) * Math.Log(m + 1e-16)).Sum();
            log.TheoreticalEntropy = -p.MultiplyElements(p.Ln()).Sum();

            return log;
        }

        public FiniteStateMCLog RunAndLog(int numOfSteps, int intervalBetweenLogEntries)
        {
            var numOfLogEntries = numOfSteps/intervalBetweenLogEntries;

            var log = new FiniteStateMCLog(numOfLogEntries);


            Stream trajectoriesStream = null;

            if (OutputStream != null)
                trajectoriesStream = OutputStream;
            else if (OutputFileName != null)
                trajectoriesStream = new FileStream(OutputFileName, FileMode.Create, FileAccess.Write, FileShare.Read);

            Action<int> writeFunction;
            BitWriter bs = null;
            var numOfBits = (int)Math.Ceiling(Math.Log(StateEnergies.Length, 2.0));

            switch (OutputStyle)
            {
                case 1:
                    writeFunction = (state => { trajectoriesStream.WriteByte((byte)(1 << state)); });
                    break;

                case 2:
                    bs = new BitWriter(trajectoriesStream);
                    writeFunction = (state => { bs.WriteBits((ulong)state, numOfBits); });
                    break;

                default:
                    writeFunction = (state => { trajectoriesStream.WriteByte((byte)state); });
                    break;

            }

            for (int stepIndex = 0; stepIndex < numOfSteps; stepIndex++)
            {
                var nextState = _random.Next(_stateEnergies.Length);

                var shouldAccept = (_stateEnergies[CurrentState] > _stateEnergies[nextState]) ||
                                   (_random.NextDouble() <
                                    Math.Exp((_stateEnergies[CurrentState] - _stateEnergies[nextState])/T));

                if (shouldAccept)
                    CurrentState = nextState;

                if (stepIndex%intervalBetweenLogEntries == 0)
                {
                    _stateVisits[CurrentState]++;

                    log.VisitLog[stepIndex / intervalBetweenLogEntries] = CurrentState;

                    // Write trajectories
                    if (trajectoriesStream != null)
                        writeFunction(CurrentState);
                }
            }

            // Dispose of the stream only if this funcion created it
            if (trajectoriesStream != null && OutputStream == null)
                trajectoriesStream.Dispose();

            double sum = _stateVisits.Sum();
            log.MeanOccupancy = _stateVisits.Select(v => v/sum).ToArray();

            log.Entropy = -log.MeanOccupancy.Select(m => (m + 1e-16)*Math.Log(m + 1e-16)).Sum();

            var p = StateEnergies.Divide(T).Negative().Exp();
            p = p.Divide(p.Sum());

            log.TheoreticalEntropy = -p.MultiplyElements(p.Ln()).Sum();

            return log;
        }
    }
}
