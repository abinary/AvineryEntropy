using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Avinery.Sim
{
    public class BinaryAssistant
    {
        const double VerySmallNumber = 1e-16;
        const double JustABitLessThanOne = (1.0 - VerySmallNumber);
        const double JustABitMoreThanOne = (1.0 + VerySmallNumber);

        public BinaryAssistant()
        {
            
        }

        public int LengthOfRandomBytePadding { get; set; }
        //public bool ShouldAddSignBitToScaleDecomposition { get; set; }

        //public class CoordinatesAlignmentLog
        //{
        //    public CoordinatesAlignmentLog(int numOfSets)
        //    {
        //        UpdatedCoordinateSets = new double[numOfSets][,];
        //        RotationMats = new Mat[numOfSets];
        //    }

        //    public int AlignmentIndex1 { get; set; }
        //    public int AlignmentIndex2 { get; set; }
        //    public int AlignmentIndex3 { get; set; }

        //    public double[][,] UpdatedCoordinateSets;
        //    public Mat[] RotationMats;
        //}

        //public static CoordinatesAlignmentLog AlignCoordinates3d(double[][,] coordinateSets, int index1 = 0, int index2 = -1, int index3 = -1)
        //{
        //    var result = new CoordinatesAlignmentLog(coordinateSets.Length);
        //    var numOfCoordinatesInSet = coordinateSets[0].GetLength(0);

        //    if (index2 == -1)
        //        index2 = numOfCoordinatesInSet - 1;

        //    if (index3 == -1)
        //        index3 = numOfCoordinatesInSet / 2;

        //    result.AlignmentIndex1 = index1;
        //    result.AlignmentIndex2 = index2;
        //    result.AlignmentIndex3 = index3;

        //    for (int i = 0; i < coordinateSets.Length; i++)
        //    {
        //        var alignmentResult = AlignCoordinates3d(coordinateSets[i], index1, index2, index3);
        //        result.UpdatedCoordinateSets[i] = alignmentResult.Item1;
        //        result.RotationMats[i] = alignmentResult.Item2;
        //    }

        //    return result;
        //}

        //public static double[,,] AlignCoordinates3d(double[,,] coordinateSets, int index1 = 0, int index2 = -1, int index3 = -1)
        //{
        //    coordinateSets = (double[,,])coordinateSets.Clone();
        //    var numOfCoordinatesInSet = coordinateSets.GetLength(1);

        //    if (index2 == -1)
        //        index2 = numOfCoordinatesInSet - 1;

        //    if (index3 == -1)
        //        index3 = numOfCoordinatesInSet / 2;

        //    AlignCoordinates3dInPlace(coordinateSets, index1, index2, index3);
        //    return coordinateSets;
        //}

        static double Phi(double[,] coordinates, int index)
        {
            return Math.Atan2(coordinates[index, 1], coordinates[index, 0]);
        }

        static double Theta(double[,] coordinates, int index)
        {
            var x = coordinates[index, 0];
            var y = coordinates[index, 1];
            return Math.Atan2(coordinates[index, 2], Math.Sqrt(x*x + y*y));
        }

        static double Phi(double[,,] coordinates, int frameIndex, int index)
        {
            return Math.Atan2(coordinates[frameIndex, index, 1], coordinates[frameIndex, index, 0]);
        }

        static double Theta(double[,,] coordinates, int frameIndex, int index)
        {
            var x = coordinates[frameIndex, index, 0];
            var y = coordinates[frameIndex, index, 1];
            return Math.Atan2(coordinates[frameIndex, index, 2], Math.Sqrt(x * x + y * y));
        }

        static double[] Row(double[,] array, int index)
        {
            var length1 = array.GetLength(1);
            var result = new double[length1];

            for (int i = 0; i < length1; i++)
                result[i] = array[index, i];

            return result;
        }

        static double[] Col(double[,] array, int index)
        {
            var length0 = array.GetLength(0);
            var result = new double[length0];

            for (int i = 0; i < length0; i++)
                result[i] = array[i, index];

            return result;
        }

        public static double[][,] ToSteps(double[][,] cartesianCoordinates)
        {
            var result = new double[cartesianCoordinates.Length][,];
            for (int i = 0; i < cartesianCoordinates.Length; i++)
                result[i] = ToSteps(cartesianCoordinates[i]);

            return result;
        }

        public static double[][,] ToStepsInSphericalCoordinates(double[][,] cartesianCoordinates, bool shouldPutSinTheta = false)
        {
            var result = new double[cartesianCoordinates.Length][,];
            for (int i = 0; i < cartesianCoordinates.Length; i++)
                result[i] = ToStepsInSphericalCoordinates(cartesianCoordinates[i], shouldPutSinTheta);

            return result;
        }

        public static double[][,] ToLocalSphericalCoordinates(double[][,] cartesianCoordinates, bool shouldPutSinTheta = false)
        {
            var result = new double[cartesianCoordinates.Length][,];
            for (int i = 0; i < cartesianCoordinates.Length; i++)
                result[i] = ToLocalSphericalCoordinates(cartesianCoordinates[i], shouldPutSinTheta);

            return result;
        }

        public static double[][,] ToLocalCoordinates(double[][,] cartesianCoordinates)
        {
            var result = new double[cartesianCoordinates.Length][,];
            for (int i = 0; i < cartesianCoordinates.Length; i++)
                result[i] = ToLocalCoordinates(cartesianCoordinates[i]);

            return result;
        }

        public static double[,] LocalToNormalCoordinates(double[,] localCoordinates, double[] origin = null)
        {
            var length0 = localCoordinates.GetLength(0);

            // each 3 points define two steps, and therefore a plane
            // the plane's normal is the cross between the two steps
            // another orthogonal axis can be extracted with the cross of that normal and the 2nd step
            // the 3 axes can be used as a basis for the next (3rd) step

            // for a tethered chain, the first step and the tethering direction (usually x) also define
            // a plane, which can be preset to be the xy plane, thus removing the rotation degree of freedom
            // I'm going to assume that the chain was already aligned to remove rotational degrees of freedom
            // so that the first coordinate is just (0, 0, 0), the 1st step is added as-is,
            // and the 2nd step is defined locally in terms of the first step, and the step/x plane etc.

            // TODO: Handle cross products with near-zero result (use previous direction?..., choose random direction?...)

            var stepVectors = new double[length0, 3];
            var axis1 = new double[length0, 3];
            var axis2 = new double[length0 - 1, 3];
            var axis3 = new double[axis1.GetLength(0) - 1, 3];
            var normalCoordinates = new double[localCoordinates.GetLength(0), localCoordinates.GetLength(1)];

            // first coordinate remains the same

            if (origin != null)
            {
                normalCoordinates[0, 0] = origin[0];
                normalCoordinates[0, 1] = origin[1];
                normalCoordinates[0, 2] = origin[2];
            }


            // 2nd local coordinate is just the step vector
            normalCoordinates[1, 0] = localCoordinates[0, 0] + localCoordinates[1, 0];
            normalCoordinates[1, 1] = localCoordinates[0, 1] + localCoordinates[1, 1];
            normalCoordinates[1, 2] = localCoordinates[0, 2] + localCoordinates[1, 2];

            stepVectors[1, 0] = localCoordinates[1, 0];
            stepVectors[1, 1] = localCoordinates[1, 1];
            stepVectors[1, 2] = localCoordinates[1, 2];

            axis1[0, 0] = 1.0; // Put (1 0 0) in first step vector (other components are already 0)

            for (var i = 0; i < length0 - 2; i++)
            {
                var invNorm = 1.0 /
                              Math.Sqrt(stepVectors[1, 0] * stepVectors[1, 0] + stepVectors[1, 1] * stepVectors[1, 1] +
                                        stepVectors[1, 2] * stepVectors[1, 2]);
                axis1[i + 1, 0] = invNorm * stepVectors[i + 1, 0];
                axis1[i + 1, 1] = invNorm * stepVectors[i + 1, 1];
                axis1[i + 1, 2] = invNorm * stepVectors[i + 1, 2];

                axis2[i, 0] = axis1[i + 1, 1] * axis1[i, 2] - axis1[i + 1, 2] * axis1[i, 1];
                axis2[i, 1] = axis1[i + 1, 2] * axis1[i, 0] - axis1[i + 1, 0] * axis1[i, 2];
                axis2[i, 2] = axis1[i + 1, 0] * axis1[i, 1] - axis1[i + 1, 1] * axis1[i, 0];

                invNorm = 1.0 / Math.Sqrt(axis2[i, 0] * axis2[i, 0] + axis2[i, 1] * axis2[i, 1] + axis2[i, 2] * axis2[i, 2]);

                axis2[i, 0] = invNorm * axis2[i, 0];
                axis2[i, 1] = invNorm * axis2[i, 1];
                axis2[i, 2] = invNorm * axis2[i, 2];

                axis3[i, 0] = axis2[i, 1] * axis1[i, 2] - axis2[i, 2] * axis1[i, 1];
                axis3[i, 1] = axis2[i, 2] * axis1[i, 0] - axis2[i, 0] * axis1[i, 2];
                axis3[i, 2] = axis2[i, 0] * axis1[i, 1] - axis2[i, 1] * axis1[i, 0];

                stepVectors[i + 2, 0] = localCoordinates[i + 2, 0] * axis1[i, 0] + localCoordinates[i + 2, 1] * axis2[i, 0] + localCoordinates[i + 2, 2] * axis3[i, 0];
                stepVectors[i + 2, 1] = localCoordinates[i + 2, 0] * axis1[i, 1] + localCoordinates[i + 2, 1] * axis2[i, 1] + localCoordinates[i + 2, 2] * axis3[i, 1];
                stepVectors[i + 2, 2] = localCoordinates[i + 2, 0] * axis1[i, 2] + localCoordinates[i + 2, 1] * axis2[i, 2] + localCoordinates[i + 2, 2] * axis3[i, 2];

                normalCoordinates[i + 2, 0] = normalCoordinates[i + 1, 0] + stepVectors[i + 2, 0];
                normalCoordinates[i + 2, 1] = normalCoordinates[i + 1, 1] + stepVectors[i + 2, 1];
                normalCoordinates[i + 2, 2] = normalCoordinates[i + 1, 2] + stepVectors[i + 2, 2];
            }

            return normalCoordinates;
            //return new double[][,] { updatedCoordinates , axis1, axis2, axis3 };
        }

        public static double[,] ToSteps(double[,] cartesianCoordinates)
        {
            var length0 = cartesianCoordinates.GetLength(0);
            var stepVectors = new double[length0, 3];

            // populate normalized step vectors
            for (int i = 1; i < length0; i++)
            {
                stepVectors[i, 0] = cartesianCoordinates[i, 0] - cartesianCoordinates[i - 1, 0];
                stepVectors[i, 1] = cartesianCoordinates[i, 1] - cartesianCoordinates[i - 1, 1];
                stepVectors[i, 2] = cartesianCoordinates[i, 2] - cartesianCoordinates[i - 1, 2];
            }

            return stepVectors;
        }

        // theta is defined as in MATLAB
        public static double[,] ToStepsInSphericalCoordinates(double[,] cartesianCoordinates, bool shouldPutSinTheta = false)
        {
            var stepVectors = ToSteps(cartesianCoordinates);
            var length0 = stepVectors.GetLength(0);

            for (int i = 0; i < length0; i++)
            {
                var rxy = Math.Sqrt(stepVectors[i, 0] * stepVectors[i, 0] + stepVectors[i, 1] * stepVectors[i, 1]);
                var theta = Math.Atan2(stepVectors[i, 2], rxy);
                var phi = Math.Atan2(stepVectors[i, 1], stepVectors[i, 0]);

                stepVectors[i, 0] = Math.Sqrt(stepVectors[i, 0] * stepVectors[i, 0] + stepVectors[i, 1] * stepVectors[i, 1] + stepVectors[i, 2] * stepVectors[i, 2]); ;
                stepVectors[i, 1] = shouldPutSinTheta ? Math.Sin(theta) : theta;
                stepVectors[i, 2] = phi;
            }

            return stepVectors;
        }

        // theta is defined as in MATLAB
        public static double[,] ToLocalSphericalCoordinates(double[,] cartesianCoordinates, bool shouldPutSinTheta = false)
        {
            var localCoordinates = ToLocalCoordinates(cartesianCoordinates);
            var length0 = localCoordinates.GetLength(0);

            for (int i = 0; i < length0; i++)
            {
                var rxy = Math.Sqrt(localCoordinates[i, 0] * localCoordinates[i, 0] + localCoordinates[i, 1] * localCoordinates[i, 1]);
                var theta = Math.Atan2(localCoordinates[i, 2], rxy);
                var phi = Math.Atan2(localCoordinates[i, 1], localCoordinates[i, 0]);

                localCoordinates[i, 0] = Math.Sqrt(localCoordinates[i, 0] * localCoordinates[i, 0] + localCoordinates[i, 1] * localCoordinates[i, 1] + localCoordinates[i, 2] * localCoordinates[i, 2]); ;
                localCoordinates[i, 1] = shouldPutSinTheta ? Math.Sin(theta) : theta;
                localCoordinates[i, 2] = phi;
            }

            return localCoordinates;
        }

        public static double[,] ToLocalCoordinates(double[,] cartesianCoordinates)
        {
            var length0 = cartesianCoordinates.GetLength(0);

            // each 3 points define two steps, and therefore a plane
            // the plane's normal is the cross between the two steps
            // another orthogonal axis can be extracted with the cross of that normal and the 2nd step
            // the 3 axes can be used as a basis for the next (3rd) step

            // for a tethered chain, the first step and the tethering direction (usually x) also define
            // a plane, which can be preset to be the xy plane, thus removing the rotation degree of freedom
            // I'm going to assume that the chain was already aligned to remove rotational degrees of freedom
            // so that the first coordinate is just (0, 0, 0), the 1st step is added as-is,
            // and the 2nd step is defined locally in terms of the first step, and the step/x plane etc.

            // TODO: Handle cross products with near-zero result (use previous direction?..., choose random direction?...)

            var stepVectors = new double[length0, 3];
            var axis1 = new double[length0, 3];
            axis1[0, 0] = 1.0; // Put (1 0 0) in first step vector (other components are already 0)

            // populate normalized step vectors
            for (int i = 1; i < length0; i++)
            {
                stepVectors[i, 0] = cartesianCoordinates[i, 0] - cartesianCoordinates[i - 1, 0];
                stepVectors[i, 1] = cartesianCoordinates[i, 1] - cartesianCoordinates[i - 1, 1];
                stepVectors[i, 2] = cartesianCoordinates[i, 2] - cartesianCoordinates[i - 1, 2];

                var invNorm = 1.0/
                              Math.Sqrt(stepVectors[i, 0]*stepVectors[i, 0] + stepVectors[i, 1]*stepVectors[i, 1] +
                                        stepVectors[i, 2]*stepVectors[i, 2]);

                axis1[i, 0] = invNorm * stepVectors[i, 0];
                axis1[i, 1] = invNorm * stepVectors[i, 1];
                axis1[i, 2] = invNorm * stepVectors[i, 2];
            }


            // 2nd axis is the cross product of two consecutive step directions
            var axis2 = new double[length0 - 1, 3];

            for (int i = 0; i < length0 - 1; i++)
            {
                axis2[i, 0] = axis1[i + 1, 1]*axis1[i, 2] - axis1[i + 1, 2]*axis1[i, 1];
                axis2[i, 1] = axis1[i + 1, 2]*axis1[i, 0] - axis1[i + 1, 0]*axis1[i, 2];
                axis2[i, 2] = axis1[i + 1, 0]*axis1[i, 1] - axis1[i + 1, 1]*axis1[i, 0];

                var invNorm = 1.0 / Math.Sqrt(axis2[i, 0] * axis2[i, 0] + axis2[i, 1] * axis2[i, 1] + axis2[i, 2] * axis2[i, 2]);

                axis2[i, 0] = invNorm * axis2[i, 0];
                axis2[i, 1] = invNorm * axis2[i, 1];
                axis2[i, 2] = invNorm * axis2[i, 2];
            }


            // 3rd axis is 2nd cross previous step direction
            var axis3 = new double[axis1.GetLength(0) - 1, 3];

            for (int i = 0; i < length0 - 1; i++)
            {
                axis3[i, 0] = axis2[i, 1]*axis1[i, 2] - axis2[i, 2]*axis1[i, 1];
                axis3[i, 1] = axis2[i, 2]*axis1[i, 0] - axis2[i, 0]*axis1[i, 2];
                axis3[i, 2] = axis2[i, 0]*axis1[i, 1] - axis2[i, 1]*axis1[i, 0];
            }

            var updatedCoordinates = (double[,])cartesianCoordinates.Clone();

            // 1st coordinate is set to (0, 0, 0)
            updatedCoordinates[0, 0] = 0;
            updatedCoordinates[0, 1] = 0;
            updatedCoordinates[0, 2] = 0;

            // 2nd coordinate is the 1st step, without the rotation around the z axis
#if true
            var r =
                Math.Sqrt(stepVectors[1, 0]*stepVectors[1, 0] + stepVectors[1, 1]*stepVectors[1, 1]);
            var theta = Math.Atan2(stepVectors[1, 2], r);
            updatedCoordinates[1, 0] = r * Math.Cos(theta);
            updatedCoordinates[1, 1] = r * Math.Sin(theta);
            updatedCoordinates[1, 2] = 0;
#else
            updatedCoordinates[1, 0] = stepVectors[1, 0];
            updatedCoordinates[1, 1] = stepVectors[1, 1];
            updatedCoordinates[1, 2] = stepVectors[1, 2];
#endif

            // 3rd+ coordinates are determined by local coordinate system defined above

            // step at index "i" is from "i-1" to "i" coordinates
            // 3rd coordinate is pointed to by i=2 step
            for (int i = 2; i < length0; i++)
            {
                updatedCoordinates[i, 0] = stepVectors[i, 0] * axis1[i - 2, 0] + stepVectors[i, 1] * axis1[i - 2, 1] + stepVectors[i, 2] * axis1[i - 2, 2];
                updatedCoordinates[i, 1] = stepVectors[i, 0] * axis2[i - 2, 0] + stepVectors[i, 1] * axis2[i - 2, 1] + stepVectors[i, 2] * axis2[i - 2, 2];
                updatedCoordinates[i, 2] = stepVectors[i, 0] * axis3[i - 2, 0] + stepVectors[i, 1] * axis3[i - 2, 1] + stepVectors[i, 2] * axis3[i - 2, 2];
            }

            return updatedCoordinates;
            //return new double[][,] { updatedCoordinates , axis1, axis2, axis3 };
        }

        public static double[][,] ToLocalCoordinates_debug(double[,] cartesianCoordinates)
        {
            var length0 = cartesianCoordinates.GetLength(0);

            // each 3 points define two steps, and therefore a plane
            // the plane's normal is the cross between the two steps
            // another orthogonal axis can be extracted with the cross of that normal and the 2nd step
            // the 3 axes can be used as a basis for the next (3rd) step

            // for a tethered chain, the first step and the tethering direction (usually x) also define
            // a plane, which can be preset to be the xy plane, thus removing the rotation degree of freedom
            // I'm going to assume that the chain was already aligned to remove rotational degrees of freedom
            // so that the first coordinate is just (0, 0, 0), the 1st step is added as-is,
            // and the 2nd step is defined locally in terms of the first step, and the step/x plane etc.

            // TODO: Handle cross products with near-zero result (use previous direction?..., choose random direction?...)

            var stepVectors = new double[length0, 3];
            var axis1 = new double[length0, 3];
            axis1[0, 0] = 1.0; // Put (1 0 0) in first step vector (other components are already 0)

            // populate normalized step vectors
            for (int i = 1; i < length0; i++)
            {
                stepVectors[i, 0] = cartesianCoordinates[i, 0] - cartesianCoordinates[i - 1, 0];
                stepVectors[i, 1] = cartesianCoordinates[i, 1] - cartesianCoordinates[i - 1, 1];
                stepVectors[i, 2] = cartesianCoordinates[i, 2] - cartesianCoordinates[i - 1, 2];

                var invNorm = 1.0 /
                              Math.Sqrt(stepVectors[i, 0] * stepVectors[i, 0] + stepVectors[i, 1] * stepVectors[i, 1] +
                                        stepVectors[i, 2] * stepVectors[i, 2]);

                axis1[i, 0] = invNorm * stepVectors[i, 0];
                axis1[i, 1] = invNorm * stepVectors[i, 1];
                axis1[i, 2] = invNorm * stepVectors[i, 2];
            }


            // 2nd axis is the cross product of two consecutive step directions
            var axis2 = new double[length0 - 1, 3];

            for (int i = 0; i < length0 - 1; i++)
            {
                axis2[i, 0] = axis1[i + 1, 1] * axis1[i, 2] - axis1[i + 1, 2] * axis1[i, 1];
                axis2[i, 1] = axis1[i + 1, 2] * axis1[i, 0] - axis1[i + 1, 0] * axis1[i, 2];
                axis2[i, 2] = axis1[i + 1, 0] * axis1[i, 1] - axis1[i + 1, 1] * axis1[i, 0];

                var invNorm = 1.0 / Math.Sqrt(axis2[i, 0] * axis2[i, 0] + axis2[i, 1] * axis2[i, 1] + axis2[i, 2] * axis2[i, 2]);

                axis2[i, 0] = invNorm * axis2[i, 0];
                axis2[i, 1] = invNorm * axis2[i, 1];
                axis2[i, 2] = invNorm * axis2[i, 2];
            }


            // 3rd axis is 2nd cross previous step direction
            var axis3 = new double[axis1.GetLength(0) - 1, 3];

            for (int i = 0; i < length0 - 1; i++)
            {
                axis3[i, 0] = axis2[i, 1] * axis1[i, 2] - axis2[i, 2] * axis1[i, 1];
                axis3[i, 1] = axis2[i, 2] * axis1[i, 0] - axis2[i, 0] * axis1[i, 2];
                axis3[i, 2] = axis2[i, 0] * axis1[i, 1] - axis2[i, 1] * axis1[i, 0];
            }

            var updatedCoordinates = (double[,])cartesianCoordinates.Clone();

            // 1st coordinate is set to (0, 0, 0)
            updatedCoordinates[0, 0] = 0;
            updatedCoordinates[0, 1] = 0;
            updatedCoordinates[0, 2] = 0;

            // 2nd coordinate is the 1st step, without the rotation around the z axis
#if false
            var r =
                Math.Sqrt(stepVectors[1, 0]*stepVectors[1, 0] + stepVectors[1, 1]*stepVectors[1, 1]);
            var theta = Math.Atan2(stepVectors[1, 2], r);
            updatedCoordinates[1, 0] = r * Math.Cos(theta);
            updatedCoordinates[1, 1] = r * Math.Sin(theta);
            updatedCoordinates[1, 2] = 0;
#else
            updatedCoordinates[1, 0] = stepVectors[1, 0];
            updatedCoordinates[1, 1] = stepVectors[1, 1];
            updatedCoordinates[1, 2] = stepVectors[1, 2];
#endif

            // 3rd+ coordinates are determined by local coordinate system defined above

            // step at index "i" is from "i-1" to "i" coordinates
            // 3rd coordinate is pointed to by i=2 step
            for (int i = 2; i < length0; i++)
            {
                updatedCoordinates[i, 0] = stepVectors[i, 0] * axis1[i - 2, 0] + stepVectors[i, 1] * axis1[i - 2, 1] + stepVectors[i, 2] * axis1[i - 2, 2];
                updatedCoordinates[i, 1] = stepVectors[i, 0] * axis2[i - 2, 0] + stepVectors[i, 1] * axis2[i - 2, 1] + stepVectors[i, 2] * axis2[i - 2, 2];
                updatedCoordinates[i, 2] = stepVectors[i, 0] * axis3[i - 2, 0] + stepVectors[i, 1] * axis3[i - 2, 1] + stepVectors[i, 2] * axis3[i - 2, 2];
            }

            //return updatedCoordinates;
            return new double[][,] { updatedCoordinates, axis1, axis2, axis3 };
        }

        //public static Tuple2<double[,], Mat> AlignCoordinates3d(double[,] coordinates, int index1, int index2, int index3)
        //{
        //    var length0 = coordinates.GetLength(0);
        //    var updatedCoordinates = (double[,])coordinates.Clone();
        //    for (int i = 0; i < length0; i++)
        //    {
        //        updatedCoordinates[i, 0] -= coordinates[index1, 0];
        //        updatedCoordinates[i, 1] -= coordinates[index1, 1];
        //        updatedCoordinates[i, 2] -= coordinates[index1, 2];
        //    }

        //    var Rz = Mat.RotationMatrix(2, -Phi(updatedCoordinates, index2));
        //    var Ry = Mat.RotationMatrix(1, -Theta(updatedCoordinates, index2));

        //    var R = Rz * Ry;

        //    var row = Row(updatedCoordinates, index3);
        //    row = row * R;

        //    var Rx = Mat.RotationMatrix(0, -Math.Atan2(row[2], row[1]));

        //    R = R * Rx;

        //    updatedCoordinates = updatedCoordinates * R;
        //    //row = Row(updatedCoordinates, index3);

        //    return new Tuple2<double[,], Mat>(updatedCoordinates, R);
        //}

        //public static void AlignCoordinates3dInPlace(double[,,] coordinates, int index1 = 0, int index2 = -1, int index3 = -1)
        //{
        //    var numOfCoordinatesInSet = coordinates.GetLength(1);

        //    if (index2 == -1)
        //        index2 = numOfCoordinatesInSet - 1;

        //    if (index3 == -1)
        //        index3 = numOfCoordinatesInSet / 2;

        //    var numOfFrames = coordinates.GetLength(0);
        //    for (int frameIndex = 0; frameIndex < numOfFrames; frameIndex++)
        //    {
        //        var frame = coordinates.Get2d(0, frameIndex);
        //        var result = AlignCoordinates3d(frame, index1, index2, index3);
        //        coordinates.Set2d(0, frameIndex, result.Item1);
        //    }
        //}

        public static double[][,] ReadConfigurationsDouble(string filepath, int length0, int length1)
        {
            var numOfCoordinatesInConfiguration = length0 * length1;

            using (var stream = new FileStream(filepath, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                var numOfConfigurations = stream.Length/(numOfCoordinatesInConfiguration * sizeof(double));
                var result = new double[numOfConfigurations][,];

                using (var reader = new BinaryReader(stream))
                {
                    for (int confIndex = 0; confIndex < numOfConfigurations; confIndex++)
                    {
                        result[confIndex] = new double[length0, length1];

                        for (int i = 0; i < length0; i++)
                            for (int j = 0; j < length1; j++)
                                result[confIndex][i, j] = reader.ReadDouble();
                    }
                }

                return result;
            }
        }

        static void Swap<T>(T[] array, int i, int j)
        {
            var tmp = array[i];
            array[i] = array[j];
            array[j] = tmp;
        }

        static void Swap<T>(T[,] array, int i, int j, int k, int l)
        {
            var tmp = array[i, j];
            array[i, j] = array[k, l];
            array[k, l] = tmp;
        }

        public static void Shuffle<T>(T[] shuffeled, int seed)
        {
            var r = new Random(seed);

            for (int i = 0; i < shuffeled.Length; i++)
            {
                var swappedIndex2 = r.Next(shuffeled.Length - i);
                var swappedIndex1 = shuffeled.Length - i - 1;

                Swap(shuffeled, swappedIndex1, swappedIndex2);
            }
        }

        public static void ShuffleConfigurations(double[][,] shuffeled, int seed = 1)
        {
            var r = new Random(seed);

            for (int i = 0; i < shuffeled.Length; i++)
            {
                var swappedIndex2 = r.Next(shuffeled.Length - i);

                var swappedIndex1 = shuffeled.Length - i - 1;
                var tmp = shuffeled[swappedIndex1];
                shuffeled[swappedIndex1] = shuffeled[swappedIndex2];
                shuffeled[swappedIndex2] = tmp;
            }

            //return shuffeled;
        }

        public static void ShuffleEachConfiguration(double[][,] shuffeled, int seed)
        {
            var r = new Random(seed);
            for (int confIndex = 0; confIndex < shuffeled.Length; confIndex++)
            {
                Shuffle(shuffeled[confIndex], r.Next()+1);
            }
        }

        public static void Shuffle(double[,] shuffeled, int seed)
        {
            var r = new Random(seed);

            var length0 = shuffeled.GetLength(0);
            var length1 = shuffeled.GetLength(1);
            var numOfElements = length0 * length1;

            for (int i = 0; i < numOfElements; i++)
            {
                var swappedIndex2 = r.Next(numOfElements - i);
                var swappedIndex1 = numOfElements - i - 1;

                Swap(shuffeled, swappedIndex1 / length1, swappedIndex1 % length1, swappedIndex2 / length1,
                    swappedIndex2 % length1);
            }
        }

        public static void Shuffle<T>(T[,] shuffeled, int seed)
        {
            var r = new Random(seed);

            var length0 = shuffeled.GetLength(0);
            var length1 = shuffeled.GetLength(1);
            var numOfElements = length0*length1;

            for (int i = 0; i < numOfElements; i++)
            {
                var swappedIndex2 = r.Next(numOfElements - i);
                var swappedIndex1 = numOfElements - i - 1;

                Swap(shuffeled, swappedIndex1/length1, swappedIndex1%length1, swappedIndex2/length1,
                    swappedIndex2%length1);
            }
        }

        public static void Shuffle1stIndex<T>(T[,] shuffeled, int seed)
        {
            var r = new Random(seed);

            var length0 = shuffeled.GetLength(0);
            var length1 = shuffeled.GetLength(1);

            for (int i = 0; i < length0; i++)
            {
                var swappedIndex2 = r.Next(length0 - i);
                var swappedIndex1 = length0 - i - 1; // Fill new order from the end

                for (int j = 0; j < length1; j++)
                    Swap(shuffeled, swappedIndex1, j, swappedIndex2, j);
            }
        }

        public static void Shuffle2ndIndex<T>(T[,] shuffeled, int seed)
        {
            var r = new Random(seed);

            var length0 = shuffeled.GetLength(0);
            var length1 = shuffeled.GetLength(1);

            for (int j = 0; j < length1; j++)
            {
                var swappedIndex2 = r.Next(length1 - j);
                var swappedIndex1 = length1 - j - 1; // Fill new order from the end

                for (int i = 0; i < length0; i++)
                    Swap(shuffeled, i, swappedIndex1, i, swappedIndex2);
            }
        }

        public static double[][,] SubtractConfigurationsMinimum(double[][,] configurations)
        {
            var means = BinaryAssistant.CalculateMinimum(configurations);
            var result = new double[configurations.Length][,];

            var length0 = configurations[0].GetLength(0);
            var length1 = configurations[0].GetLength(1);

            for (int confIndex = 0; confIndex < configurations.Length; confIndex++)
            {
                var c = (double[,])configurations[confIndex].Clone();

                for (int i = 0; i < length0; i++)
                    for (int j = 0; j < length1; j++)
                        c[i, j] -= means[i, j];

                result[confIndex] = c;
            }

            return result;
        }

        public static double[][,] SubtractConfigurationsMean(double[][,] configurations)
        {
            var means = BinaryAssistant.CalculateMeans(configurations);
            var result = new double[configurations.Length][,];

            var length0 = configurations[0].GetLength(0);
            var length1 = configurations[0].GetLength(1);

            for (int confIndex = 0; confIndex < configurations.Length; confIndex++)
            {
                var c = (double[,])configurations[confIndex].Clone();

                for (int i = 0; i < length0; i++)
                    for (int j = 0; j < length1; j++)
                        c[i, j] -= means[i, j];

                result[confIndex] = c;
            }

            return result;
        }

        public static double[,] CalculateMeans(double[][,] data)
        {
            var length0 = data[0].GetLength(0);
            var length1 = data[0].GetLength(1);

            var result = new double[length0, length1];

            for (int confIndex = 0; confIndex < data.Length; confIndex++)
            {
                var c = data[confIndex];

                for (int i = 0; i < length0; i++)
                    for (int j = 0; j < length1; j++)
                        result[i, j] += c[i, j];
            }

            var N = (double)data.Length;
            for (int i = 0; i < length0; i++)
                for (int j = 0; j < length1; j++)
                    result[i, j] /= N;

            return result;
        }

        public static double[,] CalculateStds(double[][,] data, double[,] meanValues = null)
        {
            var length0 = data[0].GetLength(0);
            var length1 = data[0].GetLength(1);

            if (meanValues == null)
                meanValues = CalculateMeans(data);

            var result = new double[length0, length1];

            for (int confIndex = 0; confIndex < data.Length; confIndex++)
            {
                var c = data[confIndex];

                for (int i = 0; i < length0; i++)
                {
                    for (int j = 0; j < length1; j++)
                    {
                        var x = (c[i, j] - meanValues[i, j]);
                        result[i, j] += x*x;
                    }
                }
            }

            var N = (double)data.Length;
            for (int i = 0; i < length0; i++)
                for (int j = 0; j < length1; j++)
                    result[i, j] = Math.Sqrt(result[i, j] / (N - 1));

            return result;
        }

        public static double[,] CalculateMinimum(double[][,] data)
        {
            var length0 = data[0].GetLength(0);
            var length1 = data[0].GetLength(1);

            var result = new double[length0, length1];

            for (int i = 0; i < length0; i++)
                for (int j = 0; j < length1; j++)
                    result[i, j] = double.MaxValue;

            for (int confIndex = 0; confIndex < data.Length; confIndex++)
            {
                var c = data[confIndex];

                for (int i = 0; i < length0; i++)
                {
                    for (int j = 0; j < length1; j++)
                    {
                        if (result[i, j] > c[i, j])
                            result[i, j] = c[i, j];
                    }
                }
            }

            return result;
        }

        public static double[][,] CalculateMinMax(double[][,] data)
        {
            var result = new double[][,] {
                    new double[data[0].GetLength(0), data[0].GetLength(1)],
                    new double[data[0].GetLength(0), data[0].GetLength(1)] };

            var length0 = data[0].GetLength(0);
            var length1 = data[0].GetLength(1);

            for (int i = 0; i < length0; i++)
            {
                for (int j = 0; j < length1; j++)
                {
                    result[0][i, j] = double.MaxValue;
                    result[1][i, j] = double.MinValue;
                }
            }

            for (int confIndex = 0; confIndex < data.Length; confIndex++)
            {
                var c = data[confIndex];

                for (int i = 0; i < length0; i++)
                {
                    for (int j = 0; j < length1; j++)
                    {
                        if (result[0][i, j] > c[i, j])
                            result[0][i, j] = c[i, j];

                        if (result[1][i, j] < c[i, j])
                            result[1][i, j] = c[i, j];
                    }
                }
            }

            return result;
        }

        public static double[] GetMaximumAbsPerSite(double[,] linearConfigurations, double[] maximums = null)
        {
            var length0 = linearConfigurations.GetLength(0);
            var length1 = linearConfigurations.GetLength(1);

            if (maximums == null)
            {
                maximums = new double[length1];

                for (int i = 0; i < length1; i++)
                    maximums[i] = double.NegativeInfinity;
            }

            for (int i = 0; i < length0; i++)
            {
                for (int j = 0; j < length1; j++)
                {
                    var v = Math.Abs(linearConfigurations[i, j]);
                    if (maximums[j] < v)
                        maximums[j] = v;
                }
            }

            return maximums;
        }

        public static double[,] GetMaximumAbsPerSite(double[][,] configurations2d, double[,] maximums = null)
        {
            var length0 = configurations2d[0].GetLength(0);
            var length1 = configurations2d[0].GetLength(1);

            if (maximums == null)
            {
                maximums = new double[length0, length1];

                for (int i = 0; i < length0; i++)
                    for (int j = 0; j < length1; j++)
                        maximums[i, j] = double.NegativeInfinity;
            }

            for (int confIndex = 0; confIndex < configurations2d.Length; confIndex++)
            {
                var c = configurations2d[confIndex];

                for (int i = 0; i < length0; i++)
                {
                    for (int j = 0; j < length1; j++)
                    {
                        var v = Math.Abs(c[i, j]);
                        if (maximums[i, j] < v)
                            maximums[i, j] = v;
                    }
                }
            }

            return maximums;
        }

        public static double[][,] DiscardAbsValuesAboveThresholdPerSite(double[][,] configurations2d, double[,] maximumsExclusive)
        {
            var result = new List<double[,]>(configurations2d.Length);
            var length0 = configurations2d[0].GetLength(0);
            var length1 = configurations2d[0].GetLength(1);

            for (int confIndex = 0; confIndex < configurations2d.Length; confIndex++)
            {
                var shouldInclude = true;
                var c = configurations2d[confIndex];

                for (int i = 0; i < length0; i++)
                {
                    for (int j = 0; j < length1; j++)
                    {
                        var v = Math.Abs(c[i, j]);
                        if (v >= maximumsExclusive[i, j])
                        {
                            shouldInclude = false;
                            break;
                        }
                    }
                }

                if (shouldInclude)
                    result.Add(c);
            }

            return result.ToArray();
        }

        public static double[,] GetMaximumPerSite(double[][,] configurations2d, double[,] maximums = null)
        {
            var length0 = configurations2d[0].GetLength(0);
            var length1 = configurations2d[0].GetLength(1);

            if (maximums == null)
            {
                maximums = new double[length0, length1];

                for (int i = 0; i < length0; i++)
                    for (int j = 0; j < length1; j++)
                        maximums[i, j] = double.NegativeInfinity;
            }

            for (int confIndex = 0; confIndex < configurations2d.Length; confIndex++)
            {
                var c = configurations2d[confIndex];

                for (int i = 0; i < length0; i++)
                {
                    for (int j = 0; j < length1; j++)
                    {
                        if (maximums[i, j] < c[i, j])
                            maximums[i, j] = c[i, j];
                    }
                }
            }

            return maximums;
        }

        public static double[,] GetMinimumPerSite(double[][,] configurations2d, double[,] minimums = null)
        {
            var length0 = configurations2d[0].GetLength(0);
            var length1 = configurations2d[0].GetLength(1);

            if (minimums == null)
            {
                minimums = new double[length0, length1];

                for (int i = 0; i < length0; i++)
                    for (int j = 0; j < length1; j++)
                        minimums[i, j] = double.PositiveInfinity;
            }

            for (int confIndex = 0; confIndex < configurations2d.Length; confIndex++)
            {
                var c = configurations2d[confIndex];

                for (int i = 0; i < length0; i++)
                {
                    for (int j = 0; j < length1; j++)
                    {
                        if (minimums[i, j] > c[i, j])
                            minimums[i, j] = c[i, j];
                    }
                }
            }

            return minimums;
        }



        /// <summary>
        /// 
        /// </summary>
        /// <param name="writer"></param>
        /// <param name="configurations2d"></param>
        /// <param name="numOfBitsPerValue">up to 64</param>
        /// <returns></returns>
        public static void WriteBinaryIntegers(BitWriter writer, ulong[][,] configurations2d,
            int numOfBitsPerValue)
        {
            for (int confIndex = 0; confIndex < configurations2d.Length; confIndex++)
                WriteBinaryIntegers(writer, configurations2d[confIndex], numOfBitsPerValue);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="writer"></param>
        /// <param name="configurations2d"></param>
        /// <param name="numOfBitsPerValue">up to 64</param>
        /// <returns></returns>
        public static void WriteBinaryIntegers(BitWriter writer, ulong[,] configurations2d,
            int numOfBitsPerValue)
        {
            var length0 = configurations2d.GetLength(0);
            var length1 = configurations2d.GetLength(1);

            double maximalInteger = ((1 << numOfBitsPerValue) - 1);
            double maximalIntegerExclusive = (1 << numOfBitsPerValue);
            double twiceMaximalIntegerExclusive = 2.0 * maximalIntegerExclusive;

            for (int i = 0; i < length0; i++)
                for (int j = 0; j < length1; j++)
                    writer.WriteBits(configurations2d[i, j], numOfBitsPerValue);
        }

        public static ulong[][,] ReadBinaryIntegers2dConfigurations(BitReader reader, int configurationLength0, int configurationLength1, int numOfBitsPerValue)
        {
            var numOfValues = reader.Length / numOfBitsPerValue;
            var configurationLength = configurationLength0 * configurationLength1;
            var numOfConfigurations = numOfValues / configurationLength;

            var result = new ulong[numOfConfigurations][,];

            for (int i = 0; i < numOfConfigurations; i++)
            {
                var c = new ulong[configurationLength0, configurationLength1];
                result[i] = c;

                for (int j = 0; j < configurationLength0; j++)
                    for (int k = 0; k < configurationLength1; k++)
                        c[j, k] = reader.ReadBits(numOfBitsPerValue);
            }

            return result;
        }

        public static ulong[,,] ReadBinaryIntegers(BitReader reader, int configurationLength0, int configurationLength1, int numOfBitsPerValue)
        {
            var numOfValues = reader.Length / numOfBitsPerValue;
            var configurationLength = configurationLength0 * configurationLength1;
            var numOfConfigurations = numOfValues / configurationLength;

            var result = new ulong[numOfConfigurations, configurationLength0, configurationLength1];

            for (int i = 0; i < numOfConfigurations; i++)
                for (int j = 0; j < configurationLength0; j++)
                    for (int k = 0; k < configurationLength1; k++)
                        result[i, j, k] = reader.ReadBits(numOfBitsPerValue);

            return result;
        }

        public static ulong[,] ReadBinaryIntegers(BitReader reader, int configurationLength, int numOfBitsPerValue)
        {
            var numOfValues = reader.Length / numOfBitsPerValue;
            var numOfConfigurations = numOfValues / configurationLength;

            var result = new ulong[numOfConfigurations, configurationLength];

            for (int i = 0; i < numOfConfigurations; i++)
                for (int j = 0; j < configurationLength; j++)
                    result[i, j] = reader.ReadBits(numOfBitsPerValue);

            return result;
        }

        public static ulong[] ReadBinaryIntegers(BitReader reader, int numOfBitsPerValue)
        {
            var numOfValues = reader.Length / numOfBitsPerValue;
            var result = new ulong[numOfValues];

            for (int i = 0; i < result.Length; i++)
                    result[i] = reader.ReadBits(numOfBitsPerValue);

            return result;
        }

        public static ulong[,] EmbedCoarseGrainedConfigurationsInNewAlphabet(ulong[,] configurationsLinear, ulong[,] newAlphabet)
        {
            var numOfConfigurations = configurationsLinear.GetLength(0);
            var numOfCoordinates = configurationsLinear.GetLength(1);
            var reassignedConfigurations = new ulong[numOfConfigurations, numOfCoordinates];

            for (int i = 0; i < numOfConfigurations; i++)
                for (int j = 0; j < numOfCoordinates; j++)
                    reassignedConfigurations[i, j] = newAlphabet[j, configurationsLinear[i, j]];

            return reassignedConfigurations;
        }

        // Write fixed number of bits from an integer representation in the given range per coordinate
        // That is - if k bits are used, the integer would be (value/range)*(2^k - 1)
        public static ulong[][,] ConvertToFixedNumOfStates(double[][,] configurations2d,
            double[,] range, int numOfStates, bool containsNegativeValues = true)
        {
            var cgConfigurations2d = new ulong[configurations2d.Length][,];
            var array = configurations2d[0];
            var length0 = array.GetLength(0);
            var length1 = array.GetLength(1);

            if (range == null)
            {
                range = GetMaximumAbsPerSite(configurations2d);

                for (int i = 0; i < length0; i++)
                    for (int j = 0; j < length1; j++)
                        range[i, j] *= JustABitMoreThanOne;
            }

            var numberOfValueBits = (int)Math.Ceiling(Math.Log(numOfStates, 2.0));

            //double maximalInteger = ((1 << numberOfValueBits) - 1);
            double maximalIntegerExclusive = numOfStates;
            double halfMaximalIntegerExclusive = 0.5 * maximalIntegerExclusive;

            var mask = ~(UInt64.MaxValue << numberOfValueBits);

            // Write sign
            for (int stepIndex = 0; stepIndex < configurations2d.Length; stepIndex++)
            {
                var cgArray = new ulong[array.GetLength(0), array.GetLength(1)];
                cgConfigurations2d[stepIndex] = cgArray;
                array = configurations2d[stepIndex];

                for (int i = 0; i < length0; i++)
                {
                    for (int j = 0; j < length1; j++)
                    {
                        // Assuming "range" takes the (inclusive) maximal values that can appear in "array", then
                        // "abs(array[i, j]) / range[i, j]" is in the inclusive range [0, 1]
                        // In that case, the factor "JustABitLessThanOne" should be used to effectively change the range to [0,1),
                        // so that when multiplied by "maximalIntegerExclusive", the value "maximalIntegerExclusive" can never be reached.
                        // It is best to handle it in the "range" instead.

                        if (containsNegativeValues)
                            cgArray[i, j] = (ulong)((array[i, j] + range[i, j]) * (halfMaximalIntegerExclusive / range[i, j]));
                        else
                            cgArray[i, j] = (ulong)(Math.Abs(array[i, j]) / range[i, j] * maximalIntegerExclusive);
                    }
                }
            }

            return cgConfigurations2d;
        }

        // Write fixed number of bits from an integer representation in the given range per coordinate
        // That is - if k bits are used, the integer would be (value/range)*(2^k - 1)
        public static double[,] WriteBinaryFixedNumOfStates(BitWriter writer, double[][,] configurations2d,
            double[,] range, int numOfStates, bool containsNegativeValues = true, int bitOutputStyle = 0)
        {
            if (range == null)
            {
                range = GetMaximumAbsPerSite(configurations2d);
                var length0 = range.GetLength(0);
                var length1 = range.GetLength(1);

                for (int i = 0; i < length0; i++)
                    for (int j = 0; j < length1; j++)
                        range[i, j] *= JustABitMoreThanOne;
            }

            var cgConfigurations = ConvertToFixedNumOfStates(configurations2d, range, numOfStates, containsNegativeValues);

            var numberOfValueBits = (int)Math.Ceiling(Math.Log(numOfStates, 2.0));

            switch (bitOutputStyle)
            {
                case 0:
                    break;

                case 1: // Output rounded to bytes
                    numberOfValueBits = 8 * (int)Math.Ceiling(numberOfValueBits / 8.0);
                    break;

                case 2: // Output rounded to 32 bits
                    numberOfValueBits = 32 * (int)Math.Ceiling(numberOfValueBits / 32.0);
                    break;
            }

            WriteBinaryIntegers(writer, cgConfigurations, numberOfValueBits);

            return range;
        }

        public static void WriteBinarySingle(string outputFilePath, double[][,] array)
        {
            using (var stream = new FileStream(outputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read))
            {
                using (var writer = new BinaryWriter(stream))
                    BinaryAssistant.WriteBinarySingle(writer, array);
            }
        }

        public static void WriteBinarySingle(Stream stream, double[][,] array)
        {
            using (var writer = new BinaryWriter(stream))
                BinaryAssistant.WriteBinarySingle(writer, array);
        }

        public static void WriteBinarySingle(BinaryWriter writer, double[][,] array)
        {
            for (int i = 0; i < array.Length; i++)
                WriteBinarySingle(writer, array[i]);
        }

        public static void WriteBinarySingle(BinaryWriter writer, double[,] array)
        {
            var length0 = array.GetLength(0);
            var length1 = array.GetLength(1);

            for (int i = 0; i < length0; i++)
            {
                for (int j = 0; j < length1; j++)
                    writer.Write((float)array[i, j]);
            }
        }

        public static void WriteBinary(string outputFilePath, double[][,] array)
        {
            using (var stream = new FileStream(outputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read))
            {
                using (var writer = new BinaryWriter(stream))
                    BinaryAssistant.WriteBinary(writer, array);
            }
        }

        public static void WriteBinary(Stream stream, double[][,] array)
        {
            using (var writer = new BinaryWriter(stream))
                BinaryAssistant.WriteBinary(writer, array);
        }

        public static void WriteBinary(BinaryWriter writer, double[][,] array)
        {
            for (int i = 0; i < array.Length; i++)
                WriteBinary(writer, array[i]);
        }

        public static void WriteBinary(string outputFilePath, double[,] array)
        {
            using (var stream = new FileStream(outputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read))
                BinaryAssistant.WriteBinary(stream, array);
        }

        public static void WriteBinary(Stream stream, double[,] array)
        {
            using (var writer = new BinaryWriter(stream))
                BinaryAssistant.WriteBinary(writer, array);
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


    }
}
