using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Avinery.Sim
{
    public static partial class NumericArrayExtensions
    {
        public static T[][,] ToArray2d<T>(this T[,,] coordinateSets)
        {
            var length0 = coordinateSets.GetLength(0);
            var length1 = coordinateSets.GetLength(1);
            var length2 = coordinateSets.GetLength(2);

            var result = new T[length0][,];

            for (int i = 0; i < length0; i++)
            {
                var a = new T[length1, length2];

                for (int j = 0; j < length1; j++)
                    for (int k = 0; k < length2; k++)
                        a[j, k] = coordinateSets[i, j, k];

                result[i] = a;
            }

            return result;
        }

        public static ulong[,,] ToArray3d(this ulong[][,] coordinateSets)
        {
            var length0 = coordinateSets.Length;
            var length1 = coordinateSets[0].GetLength(0);
            var length2 = coordinateSets[0].GetLength(1);

            var result = new ulong[length0, length1, length2];

            for (int i = 0; i < length0; i++)
            {
                var a = coordinateSets[i];
                for (int j = 0; j < length1; j++)
                    for (int k = 0; k < length2; k++)
                        result[i, j, k] = a[j, k];
            }

            return result;
        }

        public static T[,,] ToArray3d<T>(this T[][,] coordinateSets)
        {
            var length0 = coordinateSets.Length;
            var length1 = coordinateSets[0].GetLength(0);
            var length2 = coordinateSets[0].GetLength(1);

            var result = new T[length0, length1, length2];

            for (int i = 0; i < length0; i++)
            {
                var a = coordinateSets[i];
                for (int j = 0; j < length1; j++)
                    for (int k = 0; k < length2; k++)
                        result[i, j, k] = a[j, k];
            }

            return result;
        }

        public static double[][,] ToArray2d(this double[,,] coordinateSets)
        {
            var length0 = coordinateSets.GetLength(0);
            var length1 = coordinateSets.GetLength(1);
            var length2 = coordinateSets.GetLength(2);

            var result = new double[length0][,];

            for (int i = 0; i < length0; i++)
            {
                var a = new double[length1, length2];

                for (int j = 0; j < length1; j++)
                    for (int k = 0; k < length2; k++)
                        a[j, k] = coordinateSets[i, j, k];

                result[i] = a;
            }

            return result;
        }

        public static double[,,] ToArray3d(this double[][,] coordinateSets)
        {
            var length0 = coordinateSets.Length;
            var length1 = coordinateSets[0].GetLength(0);
            var length2 = coordinateSets[0].GetLength(1);

            var result = new double[length0, length1, length2];

            for (int i = 0; i < length0; i++)
            {
                var a = coordinateSets[i];
                for (int j = 0; j < length1; j++)
                    for (int k = 0; k < length2; k++)
                        result[i, j, k] = a[j, k];
            }

            return result;
        }

        public static int[][,] ToArray2d(this int[,,] coordinateSets)
        {
            var length0 = coordinateSets.GetLength(0);
            var length1 = coordinateSets.GetLength(1);
            var length2 = coordinateSets.GetLength(2);

            var result = new int[length0][,];

            for (int i = 0; i < length0; i++)
            {
                var a = new int[length1, length2];

                for (int j = 0; j < length1; j++)
                for (int k = 0; k < length2; k++)
                    a[j, k] = coordinateSets[i, j, k];

                result[i] = a;
            }

            return result;
        }

        public static int[,,] ToArray3d(this int[][,] coordinateSets)
        {
            var length0 = coordinateSets.Length;
            var length1 = coordinateSets[0].GetLength(0);
            var length2 = coordinateSets[0].GetLength(1);

            var result = new int[length0, length1, length2];

            for (int i = 0; i < length0; i++)
            {
                var a = coordinateSets[i];
                for (int j = 0; j < length1; j++)
                for (int k = 0; k < length2; k++)
                    result[i, j, k] = a[j, k];
            }

            return result;
        }



        public static double[][,] Range(this double[][,] coordinateSets, int start, int length)
        {
            var result = new double[length][,];

            for (int i = 0; i < length; i++)
                result[i] = coordinateSets[start + i];

            return result;
        }

        public static double[,] Get2d(this double[,,] coordinateSets, int dimension, int index)
        {
            double[,] result = null;

            switch (dimension)
            {
                case 0:
                    {
                        var lenth1 = coordinateSets.GetLength(1);
                        var lenth2 = coordinateSets.GetLength(2);
                        result = new double[lenth1, lenth2];

                        for (int i = 0; i < lenth1; i++)
                            for (int j = 0; j < lenth2; j++)
                                result[i, j] = coordinateSets[index, i, j];
                    }
                    break;

                case 1:
                    {
                        var lenth1 = coordinateSets.GetLength(0);
                        var lenth2 = coordinateSets.GetLength(2);
                        result = new double[lenth1, lenth2];

                        for (int i = 0; i < lenth1; i++)
                            for (int j = 0; j < lenth2; j++)
                                result[i, j] = coordinateSets[i, index, j];
                    }
                    break;

                case 2:
                    {
                        var lenth1 = coordinateSets.GetLength(0);
                        var lenth2 = coordinateSets.GetLength(1);
                        result = new double[lenth1, lenth2];

                        for (int i = 0; i < lenth1; i++)
                            for (int j = 0; j < lenth2; j++)
                                result[i, j] = coordinateSets[i, j, index];
                    }
                    break;
            }

            return result;
        }

        public static void Set2d(this double[,,] coordinateSets, int dimension, int index, double[,] values)
        {
            switch (dimension)
            {
                case 0:
                    {
                        var lenth1 = coordinateSets.GetLength(1);
                        var lenth2 = coordinateSets.GetLength(2);
                        var result = new double[lenth1, lenth2];

                        for (int i = 0; i < lenth1; i++)
                            for (int j = 0; j < lenth2; j++)
                                coordinateSets[index, i, j] = values[i, j];
                    }
                    break;

                case 1:
                    {
                        var lenth1 = coordinateSets.GetLength(0);
                        var lenth2 = coordinateSets.GetLength(2);
                        var result = new double[lenth1, lenth2];

                        for (int i = 0; i < lenth1; i++)
                            for (int j = 0; j < lenth2; j++)
                                coordinateSets[i, index, j] = values[i, j];
                    }
                    break;

                case 2:
                    {
                        var lenth1 = coordinateSets.GetLength(0);
                        var lenth2 = coordinateSets.GetLength(1);
                        var result = new double[lenth1, lenth2];

                        for (int i = 0; i < lenth1; i++)
                            for (int j = 0; j < lenth2; j++)
                                coordinateSets[i, j, index] = values[i, j];
                    }
                    break;
            }
        }

        public static double[] ToDouble(this float[] array)
        {
            var result = new double[array.Length];

            for (int i = 0; i < array.Length; i++)
                result[i] = array[i];

            return result;
        }

        public static double[,] ToDouble(this float[,] array)
        {
            var length0 = array.GetLength(0);
            var length1 = array.GetLength(1);
            var result = new double[length0, length1];

            for (int i = 0; i < length0; i++)
                for (int j = 0; j < length1; j++)
                    result[i, j] = array[i, j];

            return result;
        }

        public static double[][,] ToDouble(this float[][,] array)
        {
            var result = new double[array.Length][,];

            for (int i = 0; i < array.Length; i++)
                result[i] = ToDouble(array[i]);

            return result;
        }

        public static double[] OperateOn(this double[,] matrix, double[] vector)
        {
            double[] result = new double[vector.Length];

            for (int matrixRow = 0; matrixRow < matrix.GetLength(0); matrixRow++)
            {
                for (int i = 0; i < vector.Length; i++)
                    result[matrixRow] += matrix[matrixRow, i] * vector[i];
            }

            return result;
        }

        public static double[] Multiply(this double[,] matrix, double[] vector)
        {
            return matrix.OperateOn(vector);
        }

        public static double[,] ToColumnMatrix(this double[] vector)
        {
            double[,] result = new double[vector.Length, 1];
            for (int i = 0; i < vector.Length; i++)
                result[i, 0] = vector[i];
            return result;
        }

        public static double[,] ToRowMatrix(this double[] vector)
        {
            double[,] result = new double[1, vector.Length];
            for (int i = 0; i < vector.Length; i++)
                result[0, i] = vector[i];
            return result;
        }

        public static double[] GetRow(this double[,] matrix, int row = 0)
        {
            int len = matrix.GetLength(1);
            double[] result = new double[len];
            for (int i = 0; i < len; i++)
                result[i] = matrix[row, i];
            return result;
        }

        public static double[] GetCol(this double[,] matrix, int col = 0)
        {
            int len = matrix.GetLength(0);
            double[] result = new double[len];
            for (int i = 0; i < len; i++)
                result[i] = matrix[i, col];
            return result;
        }

        public static double Dot(this double[] v1, double[] v2)
        {
            if (v1.Length != v2.Length)
                throw new Exception("Incompatible vector sizes!");

            double result = 0;

            for (int i = 0; i < v1.Length; i++)
                result += v1[i] * v2[i];

            return result;
        }

        public static double[] Negative(this double[] v)
        {
            var result = new double[v.Length];

            for (int i = 0; i < v.Length; i++)
                result[i] = -v[i];

            return result;
        }

        public static double[] Pow(this double[] v, double p)
        {
            var result = new double[v.Length];

            for (int i = 0; i < v.Length; i++)
                result[i] = Math.Pow(v[i], p);

            return result;
        }

        public static double[] Pow(this double v, double[] p)
        {
            var result = new double[p.Length];

            for (int i = 0; i < p.Length; i++)
                result[i] = Math.Pow(v, p[i]);

            return result;
        }

        public static double[] Pow(this double[] v, double[] p)
        {
            var result = new double[v.Length];

            for (int i = 0; i < v.Length; i++)
                result[i] = Math.Pow(v[i], p[i]);

            return result;
        }

        public static double[] Exp(this double[] v)
        {
            var result = new double[v.Length];

            for (int i = 0; i < v.Length; i++)
                result[i] = Math.Exp(v[i]);

            return result;
        }

        public static double[] Ln(this double[] v)
        {
            var result = new double[v.Length];

            for (int i = 0; i < v.Length; i++)
                result[i] = Math.Log(v[i]);

            return result;
        }

        public static double[] Log10(this double[] v)
        {
            var result = new double[v.Length];

            for (int i = 0; i < v.Length; i++)
                result[i] = Math.Log10(v[i]);

            return result;
        }

        public static double[] Log2(this double[] v)
        {
            var result = new double[v.Length];
            var invLog2 = 1.0 / Math.Log(2.0);

            for (int i = 0; i < v.Length; i++)
                result[i] = Math.Log(v[i]) * invLog2;

            return result;
        }

        public static double[] AddElements(this double[] v1, double[] v2)
        {
            var result = new double[v1.Length];

            for (int i = 0; i < v1.Length; i++)
                result[i] = v1[i] + v2[i];

            return result;
        }

        public static double[] SubtractElements(this double[] v1, double[] v2)
        {
            var result = new double[v1.Length];

            for (int i = 0; i < v1.Length; i++)
                result[i] = v1[i] - v2[i];

            return result;
        }

        public static double[] MultiplyElements(this double[] v1, double[] v2)
        {
            var result = new double[v1.Length];

            for (int i = 0; i < v1.Length; i++)
                result[i] = v1[i] * v2[i];

            return result;
        }

        public static double[,] MultiplyElements(this double[,] v1, double[,] v2)
        {
            var length0 = v1.GetLength(0);
            var length1 = v1.GetLength(1);

            var result = new double[length0, length1];

            for (int i = 0; i < length0; i++)
                for (int j = 0; j < length1; j++)
                    result[i, j] = v1[i, j] * v2[i, j];

            return result;
        }

        public static double[] DivideElements(this double[] v1, double[] v2)
        {
            var result = new double[v1.Length];

            for (int i = 0; i < v1.Length; i++)
                result[i] = v1[i] / v2[i];

            return result;
        }

        public static double[,] DivideElements(this double[,] v1, double[,] v2)
        {
            var length0 = v1.GetLength(0);
            var length1 = v1.GetLength(1);

            var result = new double[length0, length1];

            for (int i = 0; i < length0; i++)
                for (int j = 0; j < length1; j++)
                    result[i, j] = v1[i, j] / v2[i, j];

            return result;
        }

        public static double Mean(this double[] v)
        {
            return Sum(v) / v.Length;
        }

        public static int IndexOfFirstGreater(this double[] v, double value)
        {
            for (int i = 0; i < v.Length; i++)
            {
                if (v[i] > value)
                    return i;
            }

            return v.Length;
        }

        public static int IndexOfFirstLesser(this double[] v, double value)
        {
            for (int i = 0; i < v.Length; i++)
            {
                if (v[i] < value)
                    return i;
            }

            return v.Length;
        }

        public static double[] CumSum(this double[] v)
        {
            var result = new double[v.Length];
            double sum = 0;

            for (int i = 0; i < v.Length; i++)
            {
                sum += v[i];
                result[i] = sum;
            }

            return result;
        }

        public static double Sum(this double[] v)
        {
            double result = 0.0;

            for (int i = 0; i < v.Length; i++)
                result += v[i];

            return result;
        }

        public static double[] Divide(this double[] v, double s)
        {
            return v.Multiply(1.0 / s);
        }

        public static double[] Multiply(this double[] v, double s)
        {
            var result = new double[v.Length];

            for (int i = 0; i < v.Length; i++)
                result[i] = v[i] * s;

            return result;
        }

        public static double Multiply(this double[] v1, double[] v2)
        {
            return v1.Dot(v2);
        }

        public static double[,] OuterProduct(this double[] v1, double[] v2)
        {
            double[,] result = new double[v1.Length, v2.Length];

            for (int i = 0; i < v1.Length; i++)
            {
                for (int j = 0; j < v2.Length; j++)
                    result[i, j] = v1[i] * v2[j];
            }

            return result;
        }

        public static double[,] Multiply(this double[,] m1, double[,] m2)
        {
            int midCount = m1.GetLength(1);
            if (midCount != m2.GetLength(0))
                throw new Exception("Incompatible matrix sizes!");

            int rowsCount = m1.GetLength(0);
            int colsCount = m2.GetLength(1);

            double[,] result = new double[rowsCount, colsCount];
            for (int row = 0; row < rowsCount; row++)
            {
                for (int col = 0; col < colsCount; col++)
                {
                    for (int i = 0; i < midCount; i++)
                        result[row, col] = m1[row, i] * m2[i, col];
                }
            }

            return result;
        }
    }
}
