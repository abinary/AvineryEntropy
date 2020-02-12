using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Text;

namespace Avinery.Sim
{
    public class IsingSquare2dMC : MCBase
    {
        public IsingSquare2dMC()
        {
            T = 1.0;
            L = 0;
            IsPeriodic = true;

            _configuration = new int[0, 0];
        }

        private int _lastIndex;
        private int _L;
        public int L { get { return _L; } private set { _L = value; _lastIndex = value - 1; } } // Number of sites 
        public int N { get { return L*L; } }

        public double H { get; set; }
        public double J { get; set; }
        public double T { get; set; }

        private double _beta;

        public double ProbabilityForShiftMove { get; set; }
        public double ProbabilityForClusterFlipMove { get; set; }

        public int SingleFlipRejectionCount { get; private set; }
        public int SingleFlipAttemptCount { get; private set; }
        public double SingleFlipRejectionRate { get { return SingleFlipRejectionCount / (double)SingleFlipAttemptCount; } }

        public int ShiftRejectionCount { get; private set; }
        public int ShiftAttemptCount { get; private set; }
        public double ShiftRejectionRate { get { return ShiftRejectionCount / (double)ShiftAttemptCount; } }

        public bool ShouldLogTrajectories { get; set; }
        public bool ShouldSaveTrajectoriesBinary { get; set; }

        public readonly List<int> SampledSpins = new List<int>();

        public string TrajectoriesLogFileName { get; set; }

        public long ClusterFlipNumOfFlippedSites { get; private set; }

        public bool IsPeriodic { get; set; }

        private int[,] _configuration;

        public void InitializeConfigurationDown(int L)
        {
            this.L = L;
            _configuration = new int[L, L];

            for (int i = 0; i < L; i++)
                for (int j = 0; j < L; j++)
                    _configuration[i, j] = -1;
        }

        public void InitializeConfigurationUp(int L)
        {
            this.L = L;
            _configuration = new int[L, L];

            for (int i = 0; i < L; i++)
                for (int j = 0; j < L; j++)
                    _configuration[i, j] = 1;
        }

        public void InitializeConfigurationRandom(int L, double threshold = 0.5)
        {
            this.L = L;
            _configuration = new int[L, L];

            for (int i = 0; i < L; i++)
            {
                for (int j = 0; j < L; j++)
                {
                    _configuration[i, j] = (_random.NextDouble() < threshold ? -1 : 1);
                }
            }
        }

        public void InitializeConfigurationByRandomGrowth(int L, double T)
        {
            this.L = L;
            _beta = 1.0/T;
            _configuration = new int[L, L];

            _configuration[0, 0] = _random.NextDouble() < 0.5 ? -1 : 1;

            for (int size = 1; size < L; size++)
            {
                for (int i = 0; i <= size; i++)
                {
                    var up = Math.Exp(-_beta * CalculateSiteEnergy(size, i, 1)); // note the a "0" in a site doesn't contribute anything to the energy
                    var down = Math.Exp(-_beta * CalculateSiteEnergy(size, i, -1));

                    var threshold = up / (up + down);
                    _configuration[size, i] = (_random.NextDouble() <= threshold) ? 1 : -1;

                    up = Math.Exp(-_beta * CalculateSiteEnergy(i, size, 1)); // note the a "0" in a site doesn't contribute anything to the energy
                    down = Math.Exp(-_beta * CalculateSiteEnergy(i, size, -1));

                    threshold = up / (up + down);
                    _configuration[i, size] = (_random.NextDouble() <= threshold) ? 1 : -1;
                    
                    // Note: The corner [size, size] is being assigned twice. Doesn't really matter
                }
            }
        }

        public void SetConfiguration(int[,] configuration)
        {
            if (configuration.GetLength(0) != configuration.GetLength(1))
                throw new Exception("Configuration must be square!");

            L = configuration.GetLength(0);

            _configuration = new int[L, L];

            for (int i = 0; i < L; i++)
            {
                for (int j = 0; j < L; j++)
                {
                    _configuration[i, j] = Math.Sign(configuration[i, j]);
                }
            }
        }

        public int[,] GetConfiguration()
        {
            return (int[,])_configuration.Clone();
        }

        void InitializeForRun()
        {
            _beta = 1.0 / T;

            _wolffClusterIncludedSites = new int[L, L];
        }

        public Exception LastExceptionDuringRun { get; private set; }

        private long _stepNum = 0;
        public void Run(int numOfSteps)
        {
            try
            {
                InitializeForRun();

                ShiftAttemptCount = 0;
                ShiftRejectionCount = 0;

                ClusterFlipNumOfFlippedSites = 0;

                SingleFlipAttemptCount = 0;
                SingleFlipRejectionCount = 0;

                if (ProbabilityForShiftMove > 0 || ProbabilityForClusterFlipMove > 0)
                {
                    for (_stepNum = 0; _stepNum < numOfSteps; _stepNum++)
                        RunSingleStepInternal();
                }
                else
                {
                    for (_stepNum = 0; _stepNum < numOfSteps; _stepNum++)
                        DoSingleFlipMove();
                }

            }
            catch (Exception ex)
            {
                LastExceptionDuringRun = ex;
                throw ex;
            }
        }

        // For performance
        private void RunStepsBare(int numOfSteps)
        {
            if (ProbabilityForShiftMove > 0 || ProbabilityForClusterFlipMove > 0)
            {
                for (int i = 0; i < numOfSteps; i++)
                    RunSingleStepInternal();
            }
            else
            {
                for (int i = 0; i < numOfSteps; i++)
                    DoSingleFlipMove();
            }
        }

        public IsingSquare2dMCLog RunAndLog(long numOfSteps, int numOfStepsBetweenCalculations)
        {
            IsingSquare2dMCLog log = null;

            try
            {
                InitializeForRun();

                var numOfEntries = ((numOfSteps + numOfStepsBetweenCalculations - 1) / numOfStepsBetweenCalculations);
                if (numOfEntries > Int32.MaxValue)
                    throw new Exception("Number of samples is way too large!");

                log = new IsingSquare2dMCLog((int)numOfEntries);

                StreamWriter trajectoriesWriter = null;

                if (TrajectoriesLogFileName != null)
                    trajectoriesWriter = new StreamWriter(TrajectoriesLogFileName);

                ClusterFlipNumOfFlippedSites = 0;

                ShiftAttemptCount = 0;
                ShiftRejectionCount = 0;

                SingleFlipAttemptCount = 0;
                SingleFlipRejectionCount = 0;

                for (_stepNum = 0; _stepNum < numOfSteps; _stepNum++)
                {
                    RunSingleStepInternal();

                    if (_stepNum % numOfStepsBetweenCalculations == 0)
                    {
                        var logIndex = _stepNum / numOfStepsBetweenCalculations;

                        log.Mean[logIndex] = CalculateMean();
                        log.Variance[logIndex] = CalculateMeanOfSquares() - log.Mean[logIndex];

                        CalculateB(log);
                        Calculate_dn(log);

                        if (ShouldLogTrajectories)
                        {
                            log.Trajectories[logIndex] = (int[,])_configuration.Clone();

                            WriteTrajectorySnapshot(trajectoriesWriter);
                        }

                        if (SampledSpins.Count > 0)
                        {
                            if (logIndex == 0)
                                log.SampledSpins = new int[numOfEntries, SampledSpins.Count];

                            for (int i = 0; i < SampledSpins.Count; i++)
                            {
                                var x = SampledSpins[i] / L;
                                var y = SampledSpins[i] % L;

                                log.SampledSpins[logIndex, i] = _configuration[x, y];
                            }
                        }
                    }
                }

                if (trajectoriesWriter != null)
                    trajectoriesWriter.Dispose();

                log.FinalizeCalculations();
            }
            catch (Exception ex)
            {
                LastExceptionDuringRun = ex;
                throw ex;
            }

            return log;
        }

        public IsingSquare2dMCLog RunAndSampleRepeatedQuenching(int numOfSamples, double highT, int numOfStepsAtHighTemp, int numOfStepsAtDesignatedTemp)
        {
            IsingSquare2dMCLog log = null;

            try
            {
                InitializeForRun();
                var beta0 = _beta;

                //if (numOfSamples > Int32.MaxValue)
                //throw new Exception("Number of samples is way too large!");

                log = new IsingSquare2dMCLog((int)numOfSamples);

                StreamWriter trajectoriesWriter = null;

                if (TrajectoriesLogFileName != null)
                    trajectoriesWriter = new StreamWriter(TrajectoriesLogFileName);

                ClusterFlipNumOfFlippedSites = 0;

                ShiftAttemptCount = 0;
                ShiftRejectionCount = 0;

                SingleFlipAttemptCount = 0;
                SingleFlipRejectionCount = 0;

                for (int sampleIndex = 0; sampleIndex < numOfSamples; sampleIndex++)
                {
                    // switch to high temp
                    _beta = 1.0 / highT;
                    RunStepsBare(numOfStepsAtHighTemp);

                    // Go back to designated temp
                    _beta = beta0;
                    RunStepsBare(numOfStepsAtDesignatedTemp);

                    log.Mean[sampleIndex] = CalculateMean();
                    log.Variance[sampleIndex] = CalculateMeanOfSquares() - log.Mean[sampleIndex];

                    CalculateB(log);
                    Calculate_dn(log);

                    if (ShouldLogTrajectories)
                    {
                        log.Trajectories[sampleIndex] = (int[,])_configuration.Clone();

                        WriteTrajectorySnapshot(trajectoriesWriter);
                    }
                }

                if (trajectoriesWriter != null)
                    trajectoriesWriter.Dispose();

                log.FinalizeCalculations();
            }
            catch (Exception ex)
            {
                LastExceptionDuringRun = ex;
                throw ex;
            }

            return log;
        }

        public IsingSquare2dMCLog RunAndSampleRepeatedAnnealing(int numOfSamples, double highT, int numOfStepsAtEachTemp, int numTempSteps, int numOfStepsAtDesignatedTemp)
        {
            IsingSquare2dMCLog log = null;

            try
            {
                InitializeForRun();
                var beta0 = _beta;

                //if (numOfSamples > Int32.MaxValue)
                //throw new Exception("Number of samples is way too large!");

                log = new IsingSquare2dMCLog((int)numOfSamples);

                StreamWriter trajectoriesWriter = null;

                if (TrajectoriesLogFileName != null)
                    trajectoriesWriter = new StreamWriter(TrajectoriesLogFileName);

                ClusterFlipNumOfFlippedSites = 0;

                ShiftAttemptCount = 0;
                ShiftRejectionCount = 0;

                SingleFlipAttemptCount = 0;
                SingleFlipRejectionCount = 0;

                for (int sampleIndex = 0; sampleIndex < numOfSamples; sampleIndex++)
                {
                    // switch to high temp
                    for (int tempIndex = 0; tempIndex < numTempSteps; tempIndex++)
                    {
                        var t = highT + (T - highT) * ((double) tempIndex / (double) numTempSteps);
                        _beta = 1.0 / t;
                        RunStepsBare(numOfStepsAtEachTemp);
                    }

                    // Go back to designated temp
                    _beta = beta0;
                    RunStepsBare(numOfStepsAtDesignatedTemp);

                    log.Mean[sampleIndex] = CalculateMean();
                    log.Variance[sampleIndex] = CalculateMeanOfSquares() - log.Mean[sampleIndex];

                    CalculateB(log);
                    Calculate_dn(log);

                    if (ShouldLogTrajectories)
                    {
                        log.Trajectories[sampleIndex] = (int[,])_configuration.Clone();

                        WriteTrajectorySnapshot(trajectoriesWriter);
                    }
                }

                if (trajectoriesWriter != null)
                    trajectoriesWriter.Dispose();

                log.FinalizeCalculations();
            }
            catch (Exception ex)
            {
                LastExceptionDuringRun = ex;
                throw ex;
            }

            return log;
        }

        public IsingSquare2dMCLog RunAndSample(int numOfSamples, int numOfStepsBetweenCalculations)
        {
            IsingSquare2dMCLog log = null;

            try
            {
                InitializeForRun();

                //if (numOfSamples > Int32.MaxValue)
                    //throw new Exception("Number of samples is way too large!");

                log = new IsingSquare2dMCLog((int)numOfSamples);

                StreamWriter trajectoriesWriter = null;

                if (TrajectoriesLogFileName != null)
                    trajectoriesWriter = new StreamWriter(TrajectoriesLogFileName);

                ClusterFlipNumOfFlippedSites = 0;

                ShiftAttemptCount = 0;
                ShiftRejectionCount = 0;

                SingleFlipAttemptCount = 0;
                SingleFlipRejectionCount = 0;

                for (int sampleIndex = 0; sampleIndex < numOfSamples; sampleIndex++)
                {
                    RunStepsBare(numOfStepsBetweenCalculations);

                    log.Mean[sampleIndex] = CalculateMean();
                    log.Variance[sampleIndex] = CalculateMeanOfSquares() - log.Mean[sampleIndex];

                    CalculateB(log);
                    Calculate_dn(log);

                    if (ShouldLogTrajectories)
                    {
                        log.Trajectories[sampleIndex] = (int[,])_configuration.Clone();

                        WriteTrajectorySnapshot(trajectoriesWriter);
                    }
                }

                if (trajectoriesWriter != null)
                    trajectoriesWriter.Dispose();

                log.FinalizeCalculations();
            }
            catch (Exception ex)
            {
                LastExceptionDuringRun = ex;
                throw ex;
            }

            return log;
        }

        private void WriteTrajectorySnapshot(StreamWriter trajectoriesWriter)
        {
            if (trajectoriesWriter != null)
            {
                if (ShouldSaveTrajectoriesBinary)
                {
                    byte b = 0;

                    for (int row = 0; row < L; row++)
                    {
                        for (int col = 0; col < L; col++)
                        {
                            b <<= 1;

                            if (_configuration[row, col] == 1)
                                b |= 1;

                            if (col%8 == 7)
                                trajectoriesWriter.BaseStream.WriteByte(b);
                        }

                        if (L%8 != 0) // Lattice is not a multiple of 8
                        {
                            trajectoriesWriter.BaseStream.WriteByte(b);
                            b = 0;
                        }
                    }
                }
                else // Text
                {
                    for (int row = 0; row < L; row++)
                    {
                        for (int col = 0; col < L; col++)
                            trajectoriesWriter.Write("{0}", _configuration[row, col] == 1 ? '*' : '_');
                        trajectoriesWriter.WriteLine();
                    }

                    trajectoriesWriter.WriteLine();
                }
            }
        }

        public class CorrelationTestLog
        {
            public CorrelationTestLog(int numOfBlocks, int L)
            {
                NumSamplesPerBlock = new int[numOfBlocks];
                SiteMean = new double[L, L];
                SiteMeanPerBlock = new double[numOfBlocks, L, L];

                MeanSusceptibilityPerBlock = new double[numOfBlocks];
                MeanSusceptibilitySqrPerBlock = new double[numOfBlocks];
            }

            public int NumOfSamples;
            public int[] NumSamplesPerBlock;
            public double[,,] SiteMeanPerBlock;
            public double[,] SiteMean;

            public double[] MeanSusceptibilityPerBlock;
            public double[] MeanSusceptibilitySqrPerBlock;

            public double SusceptibilityAverageOverBlocks;
            public double SusceptibilityErrorOverBlocks;
            //            public double MomentVarOverBlocks;

            public double MeanMoment;
            public double MeanSusceptibility;
            public double MeanSusceptibilitySqr;
            public double MeanSusceptibilityError;
            public double TauEffective;

            public double[] OnlineAutocorrelationTaus;
        }

        public CorrelationTestLog RunCorrelationTimeTest(int numOfSteps, int numOfStepsBetweenSnapshots, int numOfBlocks = 10)
        {
            var result = new CorrelationTestLog(numOfBlocks, L);

            var length0 = _configuration.GetLength(0);
            var length1 = _configuration.GetLength(1);
            double N = length0*length1;

            var susceptibilityHistory = new double[20];
            var nextSusceptibilityHistoryPosition = 0;
            var susceptibilityAutocorrelation = new double[susceptibilityHistory.Length];
            result.OnlineAutocorrelationTaus = susceptibilityAutocorrelation;

            for (_stepNum = 0; _stepNum < numOfSteps; _stepNum++)
            {
                RunSingleStepInternal();

                if (_stepNum%numOfStepsBetweenSnapshots == 0)
                {
                    var currentBlock = result.NumOfSamples % numOfBlocks;
                    result.NumOfSamples++;
                    result.NumSamplesPerBlock[currentBlock]++;

                    int sum = 0;

                    for (int i = 0; i < length0; i++)
                    {
                        for (int j = 0; j < length1; j++)
                        {
                            sum += _configuration[i, j];
                            result.SiteMean[i, j] += _configuration[i, j];
                            result.SiteMeanPerBlock[currentBlock, i, j] += _configuration[i, j];
                        }
                    }

                    var m = (sum/N);
                    result.MeanMoment += m;
                    var s = m*m;
                    result.MeanSusceptibility += s;
                    result.MeanSusceptibilitySqr += s * s;

                    result.MeanSusceptibilityPerBlock[currentBlock] += s;
                    result.MeanSusceptibilitySqrPerBlock[currentBlock] += s * s;

                    susceptibilityHistory[nextSusceptibilityHistoryPosition++] = s;
                    nextSusceptibilityHistoryPosition %= susceptibilityHistory.Length;

                    // Accumulate autocorrelation
                    if ((_stepNum/numOfStepsBetweenSnapshots) >= susceptibilityHistory.Length)
                    {
                        for (int i = 0; i < susceptibilityHistory.Length; i++)
                        {
                            susceptibilityAutocorrelation[i] += s * susceptibilityHistory[(nextSusceptibilityHistoryPosition - i + susceptibilityHistory.Length) % susceptibilityHistory.Length];
                        }
                    }
                }
            }

            result.MeanMoment /= result.NumOfSamples;
            result.MeanSusceptibility /= result.NumOfSamples;
            result.MeanSusceptibilitySqr /= result.NumOfSamples;
            result.MeanSusceptibilityError = (result.MeanSusceptibilitySqr - result.MeanSusceptibility * result.MeanSusceptibility) / result.NumOfSamples;
            result.MeanSusceptibilityError = Math.Sqrt(result.MeanSusceptibilityError);

            double sumOverBlocks = 0;
            double sumSqrOverBlocks = 0;
            for (int currentBlock = 0; currentBlock < numOfBlocks; currentBlock++)
            {
                result.MeanSusceptibilityPerBlock[currentBlock] /= result.NumSamplesPerBlock[currentBlock];
                result.MeanSusceptibilitySqrPerBlock[currentBlock] /= result.NumSamplesPerBlock[currentBlock];

                var s = result.MeanSusceptibilityPerBlock[currentBlock];
                sumOverBlocks += s;
                sumSqrOverBlocks += s*s;
            }

            result.SusceptibilityAverageOverBlocks = sumOverBlocks / numOfBlocks;
            result.SusceptibilityErrorOverBlocks = ((sumSqrOverBlocks / numOfBlocks) - result.SusceptibilityAverageOverBlocks * result.SusceptibilityAverageOverBlocks) / numOfBlocks;
            result.SusceptibilityErrorOverBlocks = Math.Sqrt(result.SusceptibilityErrorOverBlocks);

            for (int i = 0; i < length0; i++)
            {
                for (int j = 0; j < length1; j++)
                {
                    result.SiteMean[i, j] /= result.NumOfSamples;

                    for (int currentBlock = 0; currentBlock < numOfBlocks; currentBlock++)
                        result.SiteMeanPerBlock[currentBlock, i, j] /= result.NumSamplesPerBlock[currentBlock];
                }
            }

            var x = result.MeanSusceptibilityError / result.SusceptibilityErrorOverBlocks;
            result.TauEffective = 0.5 * x * x;

            var numAutocorrelationSamples = (result.NumOfSamples - susceptibilityAutocorrelation.Length + 1);
            var c0 = (susceptibilityAutocorrelation[0] / numAutocorrelationSamples - result.MeanSusceptibility * result.MeanSusceptibility);
            var correlationTimeEstimates = new double[susceptibilityAutocorrelation.Length];

            for (int i = 1; i < susceptibilityHistory.Length; i++)
            {
                var c = (susceptibilityAutocorrelation[i] / numAutocorrelationSamples - result.MeanSusceptibility * result.MeanSusceptibility);
                c /= c0;

                correlationTimeEstimates[i] = -i/Math.Log(c);
            }

            return result;
        }

        public int[][,] RunAndReturnTrajectories(int numOfSteps, int numOfStepsBetweenSnapshots)
        {
            InitializeForRun();

            var numOfEntries = numOfSteps / numOfStepsBetweenSnapshots;
            var trajectories = new int[numOfEntries][,];

            for (_stepNum = 0; _stepNum < numOfSteps; _stepNum++)
            {
                RunSingleStepInternal();

                if (_stepNum % numOfStepsBetweenSnapshots == 0)
                    trajectories[_stepNum / numOfStepsBetweenSnapshots] = (int[,])_configuration.Clone();
            }

            return trajectories;
        }

        public void RunSingleStep()
        {
            InitializeForRun();
            RunSingleStepInternal();
        }

        void RunSingleStepInternal()
        {
            var r = _random.NextDouble();

            if (_stepNum % Int32.MaxValue == 0)
                _wolffClusterIncludedSites = new int[L, L];

            _wolffClusterIncludedSitesStamp = (int)(_stepNum + 1);

            r -= ProbabilityForShiftMove;
            if (r < 0)
            {
                DoShiftMove();
                return;
            }

            r -= ProbabilityForClusterFlipMove;
            if (r < 0)
            {
                if (!IsTriangular || J >= 0)
                    DoWolffClusterFlipMove();
                else
                    throw new Exception("No cluster move for the triangular anti-ferromagnetic model");

                return;
            }

            DoSingleFlipMove();
        }

        readonly Queue<int> _wolffClusterAddedSites = new Queue<int>();
        private double _wolffClusterPadd;
        private int[,] _wolffClusterIncludedSites;
        private int _wolffClusterIncludedSitesStamp;

        void DoWolffClusterFlipMove()
        {
            _wolffClusterPadd = 1.0 - Math.Exp(-2.0 * _beta * J);
            var x = _random.Next(L);
            var y = _random.Next(L);
            var sign = _configuration[x, y];

            _wolffClusterIncludedSites[x, y] = _wolffClusterIncludedSitesStamp; // mark as included
            _configuration[x, y] *= -1; // flip
            _wolffClusterAddedSites.Enqueue(x + y * L);

            while (_wolffClusterAddedSites.Count > 0)
            {
                ClusterFlipNumOfFlippedSites++;
                var nextSite = _wolffClusterAddedSites.Dequeue();
                var nextSiteX = nextSite % L;
                var nextSiteY = nextSite / L;

                WolffClusterConsiderAndFlipNeighbors(sign, nextSiteX, nextSiteY);
            }

        }

        //private int[][] triangularClusterMoveTriangleChoices = new[][]
        //{
        //    new int[], 
        //};

        private bool[,,] triangularAntiFerroClusterMoveBonds;
        private int[,] triangularAntiFerroClusterMoveNumSatisfiedBonds;
        List<int> triangularAntiFerroListOfSitesInClusters;

        void DoTriangularAntiFerroClusterFlipMoveWangDeSterckMelko()
        {
            // If its not the following, then implement it
            // Generalized Monte Carlo loop algorithm for two-dimensional frustrated Ising model, Wang & De Sterck, 2012

            // Choose a dual node (plaquette). There are two for every node
            var dualNodeIndex = _random.Next(2 * L * L);
            var dualNodeI = dualNodeIndex / 2 % L;
            var dualNodeJ = dualNodeIndex / (2 * L);
            var dualNodeK = dualNodeIndex % 2;

            // With respect to the square representation with one diagonal interaction, we'll define the lower plaquette 
            // (triangle) to the right of every node as 0, the upper as 1.
            // 
            // In this notation, every dual node neighbors with (i,j, 1 - k), (i + 1, j, 1 - k), (i, j - 1, 1 - k) for k = 0
            // or (i,j, 1 - k), (i - 1, j, 1 - k), (i, j + 1, 1 - k) for k = 1
            // the (x-1) expressions can be substituted with (x + 1 - 2*k) or (x + 1 - 2*(1 - k))

            // Plaquette (i, j, k) includes nodes (i, j), (i, j + 1), (i + 1, j + 1) for k = 0, and (i, j), (i + 1, j), (i + 1, j + 1) for k = 1

            var dualLatticeVisitedSites = new bool[L, L, 2];
            dualLatticeVisitedSites[dualNodeI, dualNodeJ, dualNodeK] = true;


        }

        private void WolffClusterConsiderAndFlipSingleNeighbor(int sign, int x, int y)
        {
            if (sign == _configuration[x, y] && 
                _wolffClusterIncludedSites[x, y] != _wolffClusterIncludedSitesStamp &&
                _random.NextDouble() < _wolffClusterPadd)
            {
                _wolffClusterIncludedSites[x, y] = _wolffClusterIncludedSitesStamp; // mark as included
                _configuration[x, y] *= -1; // flip
                _wolffClusterAddedSites.Enqueue(x + y * L);
            }
        }

        private void WolffClusterConsiderAndFlipNeighbors(int sign, int x, int y)
        {
            var onXEdge = (x == 0 || x == _lastIndex);
            var onYEdge = (y == 0 || y == _lastIndex);
            var isEdge = (onXEdge || onYEdge);

            if (!isEdge) // bulk
            {
                WolffClusterConsiderAndFlipSingleNeighbor(sign, x, y - 1);
                WolffClusterConsiderAndFlipSingleNeighbor(sign, x, y + 1);
                WolffClusterConsiderAndFlipSingleNeighbor(sign, x - 1, y);
                WolffClusterConsiderAndFlipSingleNeighbor(sign, x + 1, y);

                if (IsTriangular)
                {
                    WolffClusterConsiderAndFlipSingleNeighbor(sign, x - 1, y - 1);
                    WolffClusterConsiderAndFlipSingleNeighbor(sign, x + 1, y + 1);
                }
            }
            else // edge
            {
                if (IsPeriodic)
                {
                    WolffClusterConsiderAndFlipSingleNeighbor(sign, x, (y - 1 + L) % L);
                    WolffClusterConsiderAndFlipSingleNeighbor(sign, x, (y + 1) % L);
                    WolffClusterConsiderAndFlipSingleNeighbor(sign, (x - 1 + L) % L, y);
                    WolffClusterConsiderAndFlipSingleNeighbor(sign, (x + 1) % L, y);

                    if (IsTriangular)
                    {
                        WolffClusterConsiderAndFlipSingleNeighbor(sign, (x - 1 + L) % L, (y - 1 + L) % L);
                        WolffClusterConsiderAndFlipSingleNeighbor(sign, (x + 1) % L, (y + 1) % L);
                    }
                }
                else // Non-periodic handling
                {
                    if (y != 0)
                        WolffClusterConsiderAndFlipSingleNeighbor(sign, x, y - 1);

                    if (y != _lastIndex)
                        WolffClusterConsiderAndFlipSingleNeighbor(sign, x, y + 1);

                    if (x != 0)
                        WolffClusterConsiderAndFlipSingleNeighbor(sign, x - 1, y);

                    if (x != _lastIndex)
                        WolffClusterConsiderAndFlipSingleNeighbor(sign, x + 1, y);

                    if (IsTriangular)
                    {
                        if (x - 1 >= 0 && y - 1 >= 0)
                            WolffClusterConsiderAndFlipSingleNeighbor(sign, x - 1, y - 1);

                        if (x + 1 <= _lastIndex && y + 1 <= _lastIndex)
                            WolffClusterConsiderAndFlipSingleNeighbor(sign, x + 1, y + 1);
                    }
                }
            }
        }

        void DoShiftMove()
        {
            ShiftAttemptCount++;

            var x = _random.Next(L);
            var y = _random.Next(L);

            int shiftX;
            int shiftY;

            if (_random.NextDouble() < 0.5)
            {
                shiftY = 0;
                shiftX = _random.NextDouble() < 0.5 ? -1 : 1;
            }
            else
            {
                shiftX = 0;
                shiftY = _random.NextDouble() < 0.5 ? -1 : 1;
            }

            var energyBefore = CalculateSiteEnergy(x, y, _configuration[x, y]);
            energyBefore += CalculateSiteEnergy((x + shiftX + L) % L, (y + shiftY + L) % L, _configuration[(x + shiftX + L) % L, (y + shiftY + L) % L]);
            energyBefore += J*_configuration[x, y]*_configuration[(x + shiftX + L)%L, (y + shiftY + L)%L]; // Remove the contribution of the interaction between the two sites

            var energyAfter = CalculateSiteEnergy(x, y, _configuration[(x + shiftX + L) % L, (y + shiftY + L) % L]);
            energyAfter += CalculateSiteEnergy((x + shiftX + L) % L, (y + shiftY + L) % L, _configuration[x, y]);
            energyAfter += 2*J; // Remove the contribution of the interaction between each site and itself

            var energyDiff = energyAfter - energyBefore;

            var shouldAccept = (energyDiff < 0) || (_random.NextDouble() < Math.Exp(-energyDiff / T));

            if (shouldAccept)
            {
                // Swap
                var tmp = _configuration[x, y];
                _configuration[x, y] = _configuration[(x + shiftX + L) % L, (y + shiftY + L) % L];
                _configuration[(x + shiftX + L)%L, (y + shiftY + L)%L] = tmp;
            }
            else
                ShiftRejectionCount++;
        }

        void DoSingleFlipMove()
        {
            SingleFlipAttemptCount++;

            var r = _random.Next(L*L);
            var x = r / L;
            var y = r % L;
            //var x = _random.Next(L);
            //var y = _random.Next(L);

            //var energyBefore = CalculateSiteEnergy(x, y, _configuration[x, y]);
            var energyAfter = CalculateSiteEnergy(x, y, -1*_configuration[x, y]);

            //var energyDiff = energyAfter - energyBefore;
            var energyDiff = 2 * energyAfter;

            var shouldAccept = (energyDiff < 0) || (_random.NextDouble() < Math.Exp(-energyDiff*_beta));

            if (shouldAccept)
                _configuration[x, y] *= -1;
            else
                SingleFlipRejectionCount++;
        }

        public void CalculateB(IsingSquare2dMCLog log)
        {
            //CalculateB1(log);
            CalculateBk(log.B1, log.B1x, 1);
            CalculateBk(log.B2, log.B2x, 2);
            CalculateBk(log.B3, log.B3x, 3);
            CalculateBk(log.B4, log.B4x, 4);
            CalculateBk(log.B5, log.B5x, 5);
            CalculateBk(log.B6, log.B6x, 6);
            CalculateBk(log.B7, log.B7x, 7);
        }

        public void CalculateB1(IsingSquare2dMCLog log)
        {
            var lastRow = L - 1;
            var lastCol = L - (2 + 1);

            var index1 = 0;
            var index2 = 0;

            // index is the binary representation:
            // [1 . .n][1 .. n]

            for (int row = 0; row < lastRow; row++)
            {
                for (int col = 0; col < lastCol; col++)
                {
                    index1 <<= 1; // shift
                    index1 &= 0x1;
                    index1 |= (_configuration[row, col] == 1 ? 1 : 0);

                    index2 <<= 1; // shift
                    index2 &= 0x1;
                    index2 |= (_configuration[row + 1, col + 1] == 1 ? 1 : 0);

                    var index = index2 | index1 << 1;
                    log.B1[index]++;
                }
            }
        }

        public void CalculateBk(int[] Bk, int[] BkX, int k)
        {
            // Cluster width is 2k-1
            if ((2*k - 1) > L)
                return;

            var lastRow = L - 1;
            var lastCol = L - (2*k + 1);
            var firstCol = k - 1;

            var index1 = 0;
            var index2 = 0;

            // index is the binary representation:
            // [1 . .n][1 .. n]

            var bitMask = ~(-1 << k);

            for (int row = 0; row < lastRow; row++)
            {
                // Initialize the index
                for (int col = 0; col < firstCol; col++)
                {
                    index1 <<= 1; // shift
                    index1 |= (_configuration[row, col] == 1 ? 1 : 0); // insert next bit

                    // Also possible:
                    // index1 ^= (int)(((uint)_configuration[row, col]) >> 31);

                    index2 <<= 1; // shift
                    index2 |= (_configuration[row, col + k] == 1 ? 1 : 0); // insert next bit
                }

                // Scan the Bk pattern across the matrix
                for (int col = firstCol; col < lastCol; col++)
                {
                    index1 <<= 1; // shift
                    index1 &= bitMask; // zero right-most bit
                    index1 |= (_configuration[row, col] == 1 ? 1 : 0); // insert next bit

                    index2 <<= 1; // shift
                    index2 &= bitMask; // zero right-most bit
                    index2 |= (_configuration[row + 1, col + k] == 1 ? 1 : 0); // insert next bit

                    var index = index2 | (index1 << k);

                    // Count
                    Bk[index]++;

                    // Count if x
                    if (_configuration[row, col + k + 1] == 1)
                        BkX[index]++;
                }
            }
        }

        // Two-Sided Bounds on the Free Energy from Local States in Monte Carlo Simulations (A. G. Schlijper, B. Smit 1989)
        void Calculate_dn_without_edges(IsingSquare2dMCLog log)
        {
            var n = log.dn_N;
            var lastRow = L - 1;
            var lastCol = L - (2 * n);
            var firstCol = n - 1;

            var index1 = 0;
            var index2 = 0;

            // index is the binary representation:
            // [1 . .n][1 .. n]
            // The cluster Hn is considered as in the paper, only horizontal instead of vertical

            var bitMask = ~(-1 << n);

            for (int row = 0; row < lastRow; row++)
            {
                // Initialize the index
                for (int col = 0; col < firstCol; col++)
                {
                    index1 <<= 1; // shift
                    index1 |= (_configuration[row, col] == 1 ? 1 : 0); // insert next bit

                    // Also possible:
                    // index1 ^= (int)(((uint)_configuration[row, col]) >> 31);

                    index2 <<= 1; // shift
                    index2 |= (_configuration[row + 1, col + n - 1] == 1 ? 1 : 0); // insert next bit
                }

                // Scan the pattern across the matrix
                for (int col = firstCol; col < lastCol; col++)
                {
                    index1 <<= 1; // shift
                    index1 &= bitMask; // zero right-most bit
                    index1 |= (_configuration[row, col] == 1 ? 1 : 0); // insert next bit

                    index2 <<= 1; // shift
                    index2 &= bitMask; // zero right-most bit
                    index2 |= (_configuration[row + 1, col + n - 1] == 1 ? 1 : 0); // insert next bit

                    var index = index2 | (index1 << n);

                    // Count
                    log.dn_H[index]++;

                    // Count
                    index = (index2 >> 1) | (index1 << n-1);
                    log.dn_H_excluding_O[index]++;
                }
            }
        }

        void Calculate_dn(IsingSquare2dMCLog log)
        {
            var n = log.dn_N;
            var lastNoneEdgeRow = L - 1;
            var lastNonEdgeCol = L - (2 * n);
            var firstNonEdgeCol = n - 1;

            if (lastNonEdgeCol < 0)
                return;

            var index1 = 0;
            var index2 = 0;

            // index is the binary representation:
            // [1 . .n][1 .. n]

            var bitMask = ~(-1 << n);
            var bitMask2 = ~(-1 << (n-1));

            #region Bulk rows

            for (int row = 0; row < lastNonEdgeCol; row++)
            {
                // Initialize the index
                for (int col = -n + 1; col < 0; col++)
                {
                    index1 <<= 1; // shift
                    index1 |= (_configuration[row, (col + L)%L] == 1 ? 1 : 0); // insert next bit

                    // Also possible:
                    // index1 ^= (int)(((uint)_configuration[row, col]) >> 31);

                    index2 <<= 1; // shift
                    index2 |= (_configuration[row + 1, col + n - 1] == 1 ? 1 : 0); // insert next bit
                }

                // Scan the pattern across the matrix
                for (int col = 0; col < lastNonEdgeCol; col++)
                {
                    index1 <<= 1; // shift
                    index1 &= bitMask; // zero right-most bit
                    index1 |= (_configuration[row, col] == 1 ? 1 : 0); // insert next bit

                    index2 <<= 1; // shift
                    index2 &= bitMask; // zero right-most bit
                    index2 |= (_configuration[row + 1, col + n - 1] == 1 ? 1 : 0); // insert next bit

                    var index = index2 | (index1 << n);

                    // Count
                    log.dn_H[index]++;

                    // Count
                    index = (index2 & bitMask2) | (index1 << n - 1);
                    log.dn_H_excluding_O[index]++;
                }

                // Handle edge columns
                for (int col = lastNonEdgeCol; col < L; col++)
                {
                    index1 <<= 1; // shift
                    index1 &= bitMask; // zero right-most bit
                    index1 |= (_configuration[row, col] == 1 ? 1 : 0); // insert next bit

                    index2 <<= 1; // shift
                    index2 &= bitMask; // zero right-most bit
                    index2 |= (_configuration[row + 1, (col + n - 1)%L] == 1 ? 1 : 0); // insert next bit

                    var index = index2 | (index1 << n);

                    // Count
                    log.dn_H[index]++;

                    // Count
                    index = (index2 & bitMask2) | (index1 << n - 1);
                    log.dn_H_excluding_O[index]++;
                }
            }

            #endregion Bulk rows

            #region Last (edge) row

            {
                // Initialize the index
                for (int col = -n + 1; col < 0; col++)
                {
                    index1 <<= 1; // shift
                    index1 |= (_configuration[L-1, (col + L) % L] == 1 ? 1 : 0); // insert next bit

                    // Also possible:
                    // index1 ^= (int)(((uint)_configuration[row, col]) >> 31);

                    index2 <<= 1; // shift
                    index2 |= (_configuration[0, col + n - 1] == 1 ? 1 : 0); // insert next bit
                }

                // Scan the pattern across the matrix
                for (int col = 0; col < lastNonEdgeCol; col++)
                {
                    index1 <<= 1; // shift
                    index1 &= bitMask; // zero right-most bit
                    index1 |= (_configuration[L-1, col] == 1 ? 1 : 0); // insert next bit

                    index2 <<= 1; // shift
                    index2 &= bitMask; // zero right-most bit
                    index2 |= (_configuration[0, col + n - 1] == 1 ? 1 : 0); // insert next bit

                    var index = index2 | (index1 << n);

                    // Count
                    log.dn_H[index]++;

                    // Count
                    index = (index2 & bitMask2) | (index1 << n - 1);
                    log.dn_H_excluding_O[index]++;
                }

                // Handle edge columns
                for (int col = lastNonEdgeCol; col < L; col++)
                {
                    index1 <<= 1; // shift
                    index1 &= bitMask; // zero right-most bit
                    index1 |= (_configuration[L-1, col] == 1 ? 1 : 0); // insert next bit

                    index2 <<= 1; // shift
                    index2 &= bitMask; // zero right-most bit
                    index2 |= (_configuration[0, (col + n - 1) % L] == 1 ? 1 : 0); // insert next bit

                    var index = index2 | (index1 << n);

                    // Count
                    log.dn_H[index]++;

                    // Count
                    index = (index2 & bitMask2) | (index1 << n - 1);
                    log.dn_H_excluding_O[index]++;
                }
            }
            #endregion Last (edge) row
        }

        public double CalculateMean()
        {
            double sum = 0.0;

            for (int i = 0; i < L; i++)
                for (int j = 0; j < L; j++)
                    sum += _configuration[i, j];

            return sum / (double) N;
        }

        public double CalculateMeanOfSquares()
        {
            double sum = 0.0;

            for (int i = 0; i < L; i++)
                for (int j = 0; j < L; j++)
                    sum += _configuration[i, j] * _configuration[i, j];

            return sum / (double)N;
        }

        public double CalculateConfigurationEnergy()
        {
            double energy = 0.0;
            for (int i = 0; i < L; i++)
                for (int j = 0; j < L; j++)
                    energy += CalculateSiteEnergy(i, j, _configuration[i, j]);

            return energy;
        }

        // Also consider top-left & bottom-right sites as interacting neighbors
        // This is makes it effectively a triangular lattice
        public bool IsTriangular { get; set; }

        double CalculateSiteEnergy(int x, int y, int siteSign)
        {
            var onXEdge = (x == 0 || x == _lastIndex);
            var onYEdge = (y == 0 || y == _lastIndex);
            var isEdge = (onXEdge || onYEdge);
            var j = - siteSign;

            int energy = 0;

            if (!isEdge) // bulk
            {
                energy += j*_configuration[x, y - 1];
                energy += j*_configuration[x, y + 1];
                energy += j*_configuration[x - 1, y];
                energy += j*_configuration[x + 1, y];

                if (IsTriangular)
                {
                    energy += j * _configuration[x - 1, y - 1];
                    energy += j * _configuration[x + 1, y + 1];
                }
            }
            else // edge
            {
                if (IsPeriodic)
                {
                    energy += j * _configuration[x, (y - 1 + L) % L];
                    energy += j * _configuration[x, (y + 1) % L];
                    energy += j * _configuration[(x - 1 + L) % L, y];
                    energy += j * _configuration[(x + 1) % L, y];

                    if (IsTriangular)
                    {
                        energy += j * _configuration[(x - 1 + L) % L, (y - 1 + L) % L];
                        energy += j * _configuration[(x + 1) % L, (y + 1) % L];
                    }
                }
                else // Non-periodic handling
                {
                    if (y != 0)
                        energy += j * _configuration[x, y - 1];

                    if (y != _lastIndex)
                        energy += j * _configuration[x, y + 1];
                    
                    if (x != 0)
                        energy += j * _configuration[x - 1, y];

                    if (x != _lastIndex)
                        energy += j * _configuration[x + 1, y];

                    if (IsTriangular)
                    {
                        if (x-1 >= 0 && y-1 >= 0)
                            energy += j * _configuration[x - 1, y - 1];

                        if (x+1 <= _lastIndex && y+1 <= _lastIndex)
                            energy += j * _configuration[x + 1, y + 1];
                    }
                }
            }

            return -H * siteSign + J * energy;
        }
    }
}
