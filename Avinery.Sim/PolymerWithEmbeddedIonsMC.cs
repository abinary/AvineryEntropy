using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;

namespace Avinery.Sim
{
    public class PolymerWithEmbeddedIonsMC : MCBase
    {
        public double R { get; private set; }
        public int N { get; set; }

        public int MinimalMoveInterval { get; set; }
        public int MaximalMoveInterval { get; set; }

        public bool IsEdgeMoveLimitedToEnd { get; set; }
        public double ProbabilityOfMovingEdges { get; set; }
        public double ProbabilityOfRotationMove { get; set; }

        public double RotationMoveMaxAngle { get; set; }

        public double MaximalExtension => (N - 1);

        private Vec2d[] _configuration;

        public Vec2d[] Configuration
        {
            get
            {
                var result = new Vec2d[_configuration.Length];
                for (int i = 0; i < result.Length; i++)
                    result[i] = _configuration[i];
                return result;
            }
        }

        public Vec2d[] AttemptedConfiguration
        {
            get
            {
                var result = new Vec2d[_attemptedConfiguration.Length];
                for (int i = 0; i < result.Length; i++)
                    result[i] = _attemptedConfiguration[i];
                return result;
            }
        }

        public double[,] ConfigurationMatrix { get { return _configuration.ToDoubleArray(); } }
        public double[,] AttemptedConfigurationMatrix { get { return _attemptedConfiguration.ToDoubleArray(); } }

        public PolymerWithEmbeddedIonsMC()
        {
            ElectrostaticInteractionDecayKappa = 0;
            ElectrostaticInteractionPrefactor = 0;
            LenardJonesRadius = 1;
            LenardJonesPrefactor = 1;
            LenardJonesPower = 6;

            RotationMoveMaxAngle = Math.PI * 0.5;

            MaxDeltaAngleEnergyPenaltySlope = 1.0 / (Math.PI / 180);
        }

        public double BondAngleHarmonicConst { get; set; }

        public double ElectrostaticInteractionDecayKappa { get; set; }
        public double ElectrostaticInteractionPrefactor { get; set; }
        public double LenardJonesRadius { get; set; }
        public double LenardJonesPrefactor { get; set; }
        public double LenardJonesPower { get; set; }

        private readonly List<Tuple<int, double>> _ionList = new List<Tuple<int, double>>();

        public void AddIonSite(int index, double charge)
        {
            _ionList.Add(Tuple.Create(index, charge));
        }

        private int[] _ionIndexes;
        private double[] _ionCharges;

        // Should be called before running a bunch of steps
        void PrepareIonInteractionData()
        {
            // Sort by index first
            _ionList.Sort((t1, t2) => t1.Item1.CompareTo(t2.Item1));

            _ionIndexes = _ionList.Select(x => x.Item1).ToArray();
            _ionCharges = _ionList.Select(x => x.Item2).ToArray();
        }

        bool AreIonSitesInRange(int i1, int i2)
        {
            if (_ionIndexes.Length == 0)
                return false;

            // Make sure the queried indexes are sorted
            if (i2 < i1)
            {
                var t = i1;
                i1 = i2;
                i2 = t;
            }

            if (i2 < _ionIndexes[0])
                return false;

            if (i1 > _ionIndexes[_ionIndexes.Length - 1])
                return false;

            for (int ion = 0; ion < _ionIndexes.Length; ion++)
            {
                if (_ionIndexes[ion] >= i1 && _ionIndexes[ion] <= i2)
                    return true;
                else if (_ionIndexes[ion] > i2)
                    return false;
            }

            return false;
        }

        private double _energyDiff;
        public double CurrentEnergyChange { get { return _energyDiff; } }

        private double _currentEnergy;
        public double CurrentEnergy { get { return _currentEnergy; } }

        public double MaxDeltaAngle { get; set; }
        public double MaxDeltaAngleEnergyPenalty { get; set; }
        public double MaxDeltaAngleEnergyPenaltySlope { get; set; }

        double CalculateConfigurationEnergy(Vec2d[] configuration, int start = -1, int end = -1)
        {
            double energy = 0;

            if (start < 0 || AreIonSitesInRange(start, end))
            {
                for (int i = 0; i < _ionIndexes.Length; i++)
                {
                    for (int j = 0; j < i; j++)
                    {
                        // Calculate distance between sites
                        var d = (configuration[_ionIndexes[i]] - configuration[_ionIndexes[j]]).Radius;

                        energy += CalculateLenardJonesInteraction(d);
                        energy += CalculateElectrostaticInteraction(d, _ionCharges[i], _ionCharges[j]);
                    }
                }
            }

            if (MaxDeltaAngleEnergyPenalty > 0)
            {
                if (start < 0)
                    start = 0;
                if (end < 0)
                    end = configuration.Length - 1;

                // expand the range backwards to include the first dihedral angle
                if (start > 0)
                    start--;

                for (int i = start + 2; i <= end; i++)
                {
                    var a = (configuration[i] - configuration[i - 1]);
                    var b = (configuration[i - 1] - configuration[i - 2]);

                    var d = Math.Abs((a.Angle - b.Angle + Math.PI) % (2.0 * Math.PI) - Math.PI);

                    if (d > MaxDeltaAngle)
                    {
                        energy += MaxDeltaAngleEnergyPenalty + MaxDeltaAngleEnergyPenaltySlope * (d - MaxDeltaAngle);
                    }
                }
            }

            if (BondAngleHarmonicConst > 0)
            {
                if (start < 0)
                    start = 0;
                if (end < 0)
                    end = configuration.Length - 1;

                // expand the range backwards to include the first dihedral angle
                if (start > 0)
                    start--;

                for (int i = start + 2; i <= end; i++)
                {
                    //var a = (configuration[i] - configuration[i - 1]).Angle;
                    //var b = (configuration[i - 1] - configuration[i - 2]).Angle;

                    //energy += BondAngleHarmonicConst * Sqr(a - b);

                    var a = (configuration[i] - configuration[i - 1]);
                    var b = (configuration[i - 1] - configuration[i - 2]);

                    var x = Sqr((a.Angle - b.Angle + Math.PI) % (2.0 * Math.PI) - Math.PI);
                    energy += BondAngleHarmonicConst * x;

                    //energy += BondAngleHarmonicConst * (1.0 - Sqr(a.Dot(b)));
                }
            }

            return energy;
        }

        double CalculateElectrostaticInteraction(double distance, double charge1, double charge2)
        {
            var x = charge1 * charge2 / Sqr(distance) 
                    * ElectrostaticInteractionPrefactor
                    * Math.Exp(-distance * ElectrostaticInteractionDecayKappa);
            return x;
        }

        double CalculateLenardJonesInteraction(double distance)
        {
            var x = distance / LenardJonesRadius;
            x = Math.Pow(x, -LenardJonesPower);
            return LenardJonesPrefactor * (x * x - x);
        }

        static double Sqr(double x)
        {
            return x * x;
        }

        static Vec2d[] InitializeConfiguration(double R, int N, double StepSize)
        {
            var configuration = new Vec2d[N];

            if (N % 2 == 0) // even number of monomers (odd number of steps)
            {
                var yHeight = Math.Sqrt(Sqr(StepSize) - Sqr((R - StepSize) / (N - 2)));
                var spacing = (R - StepSize) / (N - 2);
                for (int i = 0; i < (N - 1); i++)
                    configuration[i].x = i * spacing;

                configuration[N - 1].x = R;

                for (int i = 1; i < (N - 1); i += 2)
                    configuration[i].y = yHeight;
            }
            else // odd number of monomers (even number of steps)
            {
                var yHeight = Math.Sqrt(Sqr(StepSize) - Sqr(R / (N - 1)));
                var spacing = R / (N - 1);
                for (int i = 0; i < N; i++)
                    configuration[i].x = i * spacing;

                configuration[N - 1].x = R; // Just in case, avoid numerical summation error

                for (int i = 1; i < N; i += 2)
                    configuration[i].y = yHeight;
            }

            return configuration;
        }

        public void InitializeConfiguration(double[,] coordinates)
        {
            _configuration = new Vec2d[coordinates.GetLength(0)];
            N = _configuration.Length;

            for (int i = 0; i < _configuration.Length; i++)
                _configuration[i] = new Vec2d(coordinates[i, 0], coordinates[i, 1]);
        }

        public void InitializeConfiguration(double R)
        {
            this.R = R;

            // TODO: The code does not work for short chains (N <= 4)

            if (R <= (2.0))
            {
                var c1 = InitializeConfiguration(2, N/2, 1);
                var c2 = InitializeConfiguration(2, N - N/2 + 1, 1);

                var pivot = new Vec2d() {x = R*0.5, y = Math.Sqrt(Sqr(2.0) - 0.25*R*R)};
                var v = pivot.DirectionVector;

                c1 = c1.Rotate(v.Angle, c1[0]);
                c2 = c2.Rotate(-v.Angle, c2[0]).Add(pivot);

                _configuration = new Vec2d[N];
                for (int i = 0; i < N/2; i++)
                    _configuration[i] = c1[i];

                for (int i = 0; i < (N - N / 2); i++)
                    _configuration[i + N/2] = c2[i+1];
            }
            else
                _configuration = InitializeConfiguration(R, N, 1.0);
        }

        private readonly double TwoPi = Math.PI*2.0;

        private Vec2d[] _attemptedConfiguration;

        void PrepareForSimulation()
        {
            if (MaximalMoveInterval > _configuration.Length)
                MaximalMoveInterval = _configuration.Length;

            if (MinimalMoveInterval < 3)
                MinimalMoveInterval = 3;

            PrepareIonInteractionData();

            _attemptedConfiguration = (Vec2d[])Configuration.Clone();
            _currentEnergy = CalculateConfigurationEnergy(_configuration);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns>boolean - whether the move was accepted</returns>
        public bool RunSingleStep()
        {
            PrepareForSimulation();
            return RunSingleStepInternal();
        }

        public bool AlternateRotationMode { get; set; }

        /// <summary>
        /// 
        /// </summary>
        /// <returns>boolean - whether the move was accepted</returns>
        public bool RunSingleStepInternal()
        {
            var accepted = false;

            int i = 0, j = 0;

            var moveTypeRandomNum = _random.NextDouble();

            moveTypeRandomNum -= ProbabilityOfMovingEdges;
            if (moveTypeRandomNum < 0)
            {
                Debug.WriteLine("Edge move");

                if (IsEdgeMoveLimitedToEnd || _random.NextDouble() < 0.5) // Choose edge
                {
                    var newPosition = _configuration[N - 2] +
                                            (new Vec2d() { x = 1.0 }).Rotate(_random.NextDouble() * TwoPi);

                    _attemptedConfiguration[N - 1] = newPosition;
                    i = N - 2;
                    j = N - 1;
                }
                else
                {
                    var newPosition = _configuration[1] +
                                      (new Vec2d() { x = 1.0 }).Rotate(_random.NextDouble() * TwoPi);

                    _attemptedConfiguration[0] = newPosition;
                    i = 0;
                    j = 1;
                }
            }
            else if (moveTypeRandomNum >= 0 && moveTypeRandomNum < ProbabilityOfRotationMove)
            {
                Debug.WriteLine("Curvature move");

                // Randomly pick a 3-node sub-chain
                i = _random.Next(N - 1);

                var dTheta = (_random.NextDouble() - 0.5) * RotationMoveMaxAngle;
                var r = _configuration[i];

                //if (_random.NextDouble() > 0.5)
                {
                    if (AlternateRotationMode)
                    {
                        _attemptedConfiguration[i + 1] = _configuration[i + 1].Rotate(dTheta, r);

                        for (int k = i + 2; k < N; k++)
                            _attemptedConfiguration[k] = (_configuration[k] - _configuration[k - 1]).DirectionVector +
                                                         _attemptedConfiguration[k - 1];
                    }
                    else
                    {
                        for (int k = i + 1; k < N; k++)
                            _attemptedConfiguration[k] = _configuration[k].Rotate(dTheta, r);
                    }

                    j = N - 1;
                }
                //else
                //{
                //    for (int k = 0; k < i; k++)
                //        _attemptedConfiguration[k] = _configuration[k].Rotate(dTheta, r);

                //    j = i;
                //    i = 0;
                //}
            }
            else
            {
                // Pick a range, flip it

                var spaceRange = MaximalMoveInterval - MinimalMoveInterval;

                int subChainSize = MinimalMoveInterval + (spaceRange > 0 ? _random.Next(spaceRange + 1) : 0);

                i = _random.Next(N - subChainSize + 1); // smallest index is 0, largest index would be N-space
                j = i + subChainSize - 1; // assuming space >= 2, smallest index is 1, largest index is N-1

                var vectorBetweenChosenSites = _configuration[j] - _configuration[i];
                var distanceBetweenChosenSites = vectorBetweenChosenSites.Norm();

                var l = vectorBetweenChosenSites.DirectionVector;
                var n = l.OrthogonalVector;
                var c0 = _configuration[i].Dot(n);

                if (subChainSize == 3)
                {
                    Debug.WriteLine("Flipping chain of size *3*");

                    // calculate the flip differently, to numerically stabilize the step size over time

                    var x = distanceBetweenChosenSites * 0.5;

                    // When the chain is very stretched, the distance can be slightly larger than
                    // 2 *StepSize due to accumulated numeric error

                    // TODO: Handle the case where x > StepSize and x < StepSize * 1.001

                    if (x < 1.0)
                    {
                        var y = Math.Sqrt(1.0 - x * x);

                        var c = _configuration[i + 1].Dot(n);
                        _attemptedConfiguration[i + 1] = _configuration[i] + x * l - Math.Sign((c - c0)) * y * n;
                    }
                    else if (x > 1.001)
                        throw new Exception("Distance between neighbors diverged!");
                }
                else
                {
                    Debug.WriteLine("Flipping chain of size {0}", subChainSize);

                    // Flip across the vector between sites by subtracting the normal component twice

                    for (int k = i + 1; k < j; k++)
                    {
                        var c = _configuration[k].Dot(n);
                        _attemptedConfiguration[k] = _configuration[k] - (2.0 * (c - c0)) * n;
                    }
                }

                // adjust to the actual range of modified positions (range should include modified bond angles, so no "j--")
                i++;
            }

            // All the moves only affect the dihedral angles of the edges. This allows some optimization
            if (_ionIndexes.Length == 0 && (j - i) >= 4)
            {
                _energyDiff =
                    CalculateConfigurationEnergy(_attemptedConfiguration, i, i + 2) -
                    CalculateConfigurationEnergy(_configuration, i, i + 2);

                _energyDiff +=
                    CalculateConfigurationEnergy(_attemptedConfiguration, j - 1, j) -
                    CalculateConfigurationEnergy(_configuration, j - 1, j);
            }
            else
            {
                _energyDiff =
                    CalculateConfigurationEnergy(_attemptedConfiguration, i, j) -
                    CalculateConfigurationEnergy(_configuration, i, j);
            }

            accepted = (_energyDiff < 0) ||
                        (_random.NextDouble() <= Math.Exp(-_energyDiff));

            if (accepted)
                _currentEnergy += _energyDiff;

            if (accepted)
            {
                // updte configurations
                for (int k = i; k <= j; k++)
                    _configuration[k] = _attemptedConfiguration[k];
            }
            else
            {
                // Revert to current configuration for next attempt
                for (int k = i; k <= j; k++)
                    _attemptedConfiguration[k] = _configuration[k];
            }

            return accepted;
        }

        [Serializable]
        public class RunLog
        {
            public RunLog(int numOfEntries)
            {
                Configurations = new double[numOfEntries][,];
                EndToEndDistances = new double[numOfEntries];
                EndToEndX = new double[numOfEntries];
                EndToEndY = new double[numOfEntries];
                Energy = new double[numOfEntries];
                Accepted = new bool[numOfEntries];
            }

            public int NumOfMoves { get; set; }
            public double[][,] Configurations { get; set; }

            public double[] EndToEndDistances { get; set; }
            public double[] EndToEndX { get; set; }
            public double[] EndToEndY { get; set; }
            public double[] Energy { get; set; }
            public bool[] Accepted { get; set; }

            public bool ShouldWriteShuffledTimeSteps { get; set; }

            public void WriteBinary(string outputFilePath)
            {
                using (var output = new FileStream(outputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read))
                {
                    using (var writer = new BinaryWriter(output))
                        BinaryAssistant.WriteBinary(writer, Configurations);
                }
            }

            public class GridEntropyResult
            {
                public double Entropy { get; set; }
            }

            public GridEntropyResult CalculateGridEntropy()
            {
                var result = new GridEntropyResult();

                var length0 = Configurations[0].GetLength(0);
                var length1 = Configurations[0].GetLength(1);

                var minMax = new Tuple2<double, double>[length1];
                for (int j = 0; j < length1; j++)
                {
                    minMax[j].Item1 = double.MaxValue;
                    minMax[j].Item2 = double.MinValue;
                }

                for (int timeStepIndex = 0; timeStepIndex < Configurations.Length; timeStepIndex++)
                {
                    var configuration = Configurations[timeStepIndex];

                    for (int i = 0; i < length0; i++)
                    {
                        for (int j = 0; j < length1; j++)
                        {
                            if (minMax[j].Item1 > configuration[i, j])
                                minMax[j].Item1 = configuration[i, j];

                            if (minMax[j].Item2 < configuration[i, j])
                                minMax[j].Item2 = configuration[i, j];
                        }
                    }
                }

                return result;
            }
        }

        public void Burn(int numOfSteps)
        {
            PrepareForSimulation();

            for (int i = 0; i < numOfSteps; i++)
            {
                RunSingleStepInternal();

                // Once in a while re-calculate entire energy (steps sum up differences, for speed)
                if (i % N == 0)
                    _currentEnergy = CalculateConfigurationEnergy(_configuration);
            }
        }

        public bool ShouldSampleConfigurations { get; set; }
        public double OutputFixedUnit { get; set; }
        public int OutputType { get; set; }
        public string TrajectoriesLogFileName { get; set; }

        private int breakAtStep = 2327574;
        public RunLog Run(int numOfSteps, int samplingInterval = 0)
        {
            PrepareForSimulation();

            var log = new RunLog(samplingInterval > 0 ? numOfSteps / samplingInterval : 0);

            BitWriter trajectoriesWriter = null;
            FileStream trajectoriesStream = null;

            if (TrajectoriesLogFileName != null)
            {
                trajectoriesStream = new FileStream(TrajectoriesLogFileName, FileMode.Create, FileAccess.Write, FileShare.Read);
                trajectoriesWriter = new BitWriter(trajectoriesStream);
            }

            for (int i = 0; i < numOfSteps; i++)
            {
                log.NumOfMoves++;

                // for debug
                if (i == breakAtStep)
                {
                }

                var accepted = RunSingleStepInternal();

                // Once in a while re-calculate entire energy (steps sum up differences, for speed)
                if (i % N == 0)
                    _currentEnergy = CalculateConfigurationEnergy(_configuration);

                if (samplingInterval > 0 && i % samplingInterval == 0)
                {
                    if (ShouldSampleConfigurations)
                        log.Configurations[i / samplingInterval] = Configuration.ToDoubleArray();

                    log.EndToEndDistances[i/samplingInterval] = (_configuration[N - 1] - _configuration[0]).Radius;
                    log.EndToEndX[i / samplingInterval] = (_configuration[N - 1] - _configuration[0]).x;
                    log.EndToEndY[i / samplingInterval] = (_configuration[N - 1] - _configuration[0]).y;
                    log.Energy[i / samplingInterval] = _currentEnergy;
                    log.Accepted[i / samplingInterval] = accepted;

                    // Write trajectories
                    if (trajectoriesWriter != null)
                    {
                        switch (OutputType)
                        {
                            case 0:
                                WriteBinary(trajectoriesWriter, Configuration.ToDoubleArray());
                                break;

                            case 1:
                                WriteBinaryFixedOutputUnit(trajectoriesWriter, Configuration.ToDoubleArray(), OutputFixedUnit);
                                break;
                        }
                    }
                }
            }

            if (trajectoriesWriter != null)
            {
                //if (OutputType == 2)
                //    log.WriteBinaryScaleDecomposedOutput(trajectoriesWriter);

                trajectoriesWriter.Dispose();
                trajectoriesStream.Dispose();
            }

            return log;
        }


        public static void WriteBinaryFixedOutputUnit(BinaryWriter writer, double[,] array, double outputUnit)
        {
            var inverseOutputUnit = 1.0 / outputUnit;
            var length0 = array.GetLength(0);
            var length1 = array.GetLength(1);

            for (int i = 0; i < length0; i++)
            {
                for (int j = 0; j < length1; j++)
                    writer.Write((long)(array[i, j] * inverseOutputUnit));
            }
        }

        public static void WriteBinary(BinaryWriter writer, double[,] array)
        {
            var length0 = array.GetLength(0);
            var length1 = array.GetLength(1);

            for (int i = 0; i < length0; i++)
            {
                for (int j = 0; j < length1; j++)
                    writer.Write(array[i, j]);
            }
        }

        private static double VerySmallNumber = 1e-10;
        public static Vec2d[] IntersectCircles(Vec2d center1, double radius1, Vec2d center2, double radius2)
        {
            var r3 = (center2 - center1).Norm();

            if (r3 > radius1 + radius2) // too far apart to intersect
                return null;

            if (radius2 > r3 + radius1) // circle one completely inside circle two
                return null;

            if (radius1 > r3 + radius2) // circle two completely inside circle one
                return null;

            if ((r3 < VerySmallNumber) & (Math.Abs(radius1 - radius2) < VerySmallNumber)) // circles identical
                return null;

            var a0 = Math.Atan2((center2.y - center1.y), (center2.x - center1.x));

            // Law of cosines
            var a1 = Math.Acos(-(Sqr(radius2) - Sqr(radius1) - Sqr(r3)) / (2 * radius1 * r3));

            var alpha1 = a0 + a1;
            var alpha2 = a0 - a1;

            var r = new Vec2d() { x = radius1 };
            var v1 = r.Rotate(alpha1);
            var v2 = r.Rotate(alpha2);

            var solution1 = center1 + v1;
            var solution2 = center1 + v2;
            return new Vec2d[] { solution1, solution2, };
        }
    }
}
