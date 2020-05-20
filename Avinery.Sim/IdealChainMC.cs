using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace Avinery.Sim
{
    public class IdealChainMC : MCBase
    {
        public double R { get; private set; }
        public int N { get; set; }
        public double PullForce { get; set; }

        public int MinimalMoveInterval { get; set; }
        public int MaximalMoveInterval { get; set; }

        public bool IsEdgeMoveLimitedToEnd { get; set; }
        public double ProbabilityOfMovingEdges { get; set; }

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

        public double[,] ConfigurationMatrix { get { return _configuration.ToDoubleArray(); } }

        public IdealChainMC()
        {
        }

        static double Sqr(double x)
        {
            return x*x;
        }

        static Vec2d[] InitializeConfiguration(double R, int N)
        {
            var configuration = new Vec2d[N];

            if (N % 2 == 0) // even number of monomers (odd number of steps)
            {
                var yHeight = Math.Sqrt(1.0 - Sqr((R - 1.0) / (N - 2)));
                var spacing = (R - 1.0) / (N - 2);
                for (int i = 0; i < (N - 1); i++)
                    configuration[i].x = i * spacing;

                configuration[N - 1].x = R;

                for (int i = 1; i < (N - 1); i += 2)
                    configuration[i].y = yHeight;
            }
            else // odd number of monomers (even number of steps)
            {
                var yHeight = Math.Sqrt(1.0 - Sqr(R / (N - 1)));
                var spacing = R / (N - 1);
                for (int i = 0; i < N; i++)
                    configuration[i].x = i * spacing;

                configuration[N - 1].x = R; // Just in case, avoid numerical summation error

                for (int i = 1; i < N; i += 2)
                    configuration[i].y = yHeight;
            }

            return configuration;
        }

        public void InitializeConfiguration(double R)
        {
            this.R = R;

            // TODO: The code does not work for short chains (N <= 4)

            if (R <= (2.0))
            {
                var c1 = InitializeConfiguration(2.0, N/2);
                var c2 = InitializeConfiguration(2.0, N - N/2 + 1);

                var pivot = new Vec2d() {x = R*0.5, y = Math.Sqrt(4.0 - 0.25*R*R)};
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
                _configuration = InitializeConfiguration(R, N);
        }

        private readonly double TwoPi = Math.PI*2.0;

        /// <summary>
        /// 
        /// </summary>
        /// <returns>boolean - whether the move was accepted</returns>
        public void RunSingleStep()
        {
            if (ProbabilityOfMovingEdges > 0 && _random.NextDouble() < ProbabilityOfMovingEdges)
            {
                //if (IsEdgeMoveLimitedToEnd)
                //{
                //    // Move only the last site, and only on X axis

                //    // Move the edge as usual
                //    _configuration[N - 1] = _configuration[N - 2] +
                //                            (new Vec2d() { x = StepSize }).Rotate(_random.NextDouble() * TwoPi);

                //    // realign the configuration with the x axis
                //    var vectorBetweenEdges = _configuration[N - 1] - _configuration[0];
                //    var angle = vectorBetweenEdges.Angle;
                //    _configuration = _configuration.Rotate(-angle, _configuration[0]);
                //}
                //else
                {
                    if (IsEdgeMoveLimitedToEnd || _random.NextDouble() < 0.5) // Choose edge
                    {
                        var newPosition = _configuration[N - 2] +
                                                (new Vec2d() { x = 1.0 }).Rotate(_random.NextDouble() * TwoPi);

                        if (PullForce == 0)
                            _configuration[N - 1] = newPosition;
                        else
                        {
                            var currentR = _configuration[N - 1].Distance(_configuration[0]);
                            var proposedR = newPosition.Distance(_configuration[0]);

                            if (proposedR >= currentR || _random.NextDouble() < Math.Exp(PullForce * (proposedR - currentR)))
                                _configuration[N - 1] = newPosition;
                        }
                    }
                    else
                    {
                        _configuration[0] = _configuration[1] +
                                            (new Vec2d() { x = 1.0 }).Rotate(_random.NextDouble() * TwoPi);
                    }
                }

                return;
            }

            // Pick a range, flip it

            var spaceRange = MaximalMoveInterval - MinimalMoveInterval;

            int subChainSize = MinimalMoveInterval + (spaceRange > 0 ? _random.Next(spaceRange + 1) : 0);

            var i = _random.Next(N - subChainSize + 1); // smallest index is 0, largest index would be N-space
            var j = i + subChainSize - 1; // assuming space >= 2, smallest index is 1, largest index is N-1

            var vectorBetweenChosenSites = _configuration[j] - _configuration[i];
            var distanceBetweenChosenSites = vectorBetweenChosenSites.Norm();

            var l = vectorBetweenChosenSites.DirectionVector;
            var n = l.OrthogonalVector;
            var c0 = _configuration[i].Dot(n);

            if (subChainSize == 3)
            {
                // calculate the flip differently, to numerically stabilize the step size over time

                var x = distanceBetweenChosenSites * 0.5;

                // When the chain is very stretched, the distance can be slightly larger than
                // 2 *StepSize due to accumulated numeric error

                // TODO: Handle the case where x > StepSize and x < StepSize * 1.001

                if (x < 1.0)
                {
                    var y = Math.Sqrt(1.0 - x*x);

                    var c = _configuration[i + 1].Dot(n);
                    _configuration[i + 1] = _configuration[i] + x*l - Math.Sign((c - c0))*y*n;
                }
                else if (x > (1.001))
                    throw new Exception("Distance between neighbors diverged!");

                return;
            }


            // Flip across the vector between sites by subtracting the normal component twice

            for (int k = i + 1; k < j; k++)
            {
                var c = _configuration[k].Dot(n);
                _configuration[k] -= (2.0 * (c - c0)) * n;
            }
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
            }

            public int NumOfMoves { get; set; }
            public double[][,] Configurations { get; set; }

            public double[] EndToEndDistances { get; set; }
            public double[] EndToEndX { get; set; }
            public double[] EndToEndY { get; set; }

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
            if (MaximalMoveInterval > _configuration.Length)
                MaximalMoveInterval = _configuration.Length;

            if (MinimalMoveInterval < 3)
                MinimalMoveInterval = 3;

            for (int i = 0; i < numOfSteps; i++)
                RunSingleStep();
        }

        public bool ShouldSampleConfigurations { get; set; }
        public double OutputFixedUnit { get; set; }
        public int OutputType { get; set; }
        public string TrajectoriesLogFileName { get; set; }

        private int breakAtStep = 2327574;
        public RunLog Run(int numOfSteps, int samplingInterval = 0)
        {
            if (MaximalMoveInterval > _configuration.Length)
                MaximalMoveInterval = _configuration.Length;

            if (MinimalMoveInterval < 3)
                MinimalMoveInterval = 3;

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

                if (i == breakAtStep)
                {
                    
                }

                RunSingleStep();

                if (samplingInterval > 0 && i % samplingInterval == 0)
                {
                    if (ShouldSampleConfigurations)
                        log.Configurations[i / samplingInterval] = Configuration.ToDoubleArray();

                    log.EndToEndDistances[i/samplingInterval] = (_configuration[N - 1] - _configuration[0]).Radius;
                    log.EndToEndX[i / samplingInterval] = (_configuration[N - 1] - _configuration[0]).x;
                    log.EndToEndY[i / samplingInterval] = (_configuration[N - 1] - _configuration[0]).y;

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
