using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Spatial.Euclidean;
using MathNet.Spatial.Units;

namespace Avinery.Sim
{
    public static class Extensions
    {
        public static double[] ToDoubleArray(this Vector2D v)
        {
            return new double[2] { v.X, v.Y };
        }

        public static double[,] ToDoubleArray(this Vector2D[] v)
        {
            var result = new double[v.Length, 2];

            for (int i = 0; i < v.Length; i++)
            {
                result[i, 0] = v[i].X;
                result[i, 1] = v[i].Y;
            }

            return result;
        }

        public static Vector2D[] ToVector2dArray(this double[,] v)
        {
            var result = new Vector2D[v.GetLength(0)];

            for (int i = 0; i < result.Length; i++)
                result[i] = new Vector2D(v[i, 0], v[i, 1]);

            return result;
        }

        public static double[] ToDoubleArray(this Vector3D v)
        {
            return new double[3] { v.X, v.Y, v.Z };
        }

        public static double[,] ToDoubleArray(this Vector3D[] v)
        {
            var result = new double[v.Length, 3];

            for (int i = 0; i < v.Length; i++)
            {
                result[i, 0] = v[i].X;
                result[i, 1] = v[i].Y;
                result[i, 2] = v[i].Z;
            }

            return result;
        }

        public static Vector3D[] ToVector3dArray(this double[,] v)
        {
            var result = new Vector3D[v.GetLength(0)];

            if (v.GetLength(1) == 3)
            {
                for (int i = 0; i < result.Length; i++)
                    result[i] = new Vector3D(v[i, 0], v[i, 1], v[i, 2]);
            }
            else
            {
                for (int i = 0; i < result.Length; i++)
                    result[i] = new Vector3D(v[i, 0], v[i, 1], 0.0);
            }

            return result;
        }

        public static Vector3D[] Add(this Vector3D[] vectors, Vector3D v)
        {
            var result = new Vector3D[vectors.Length];

            for (int i = 0; i < vectors.Length; i++)
                result[i] = vectors[i] + v;

            return result;
        }

        public static Vector3D RotatePhi(this Vector3D v, double phi)
        {
            return v.Rotate(UnitVector3D.ZAxis, Angle.FromRadians(phi));
        }

        public static Vector3D RotatePhi(this Vector3D v, Vector3D relativeTo, double phi)
        {
            return (v - relativeTo).RotatePhi(phi) + relativeTo;
        }

        public static Vector3D Rotate(this Vector3D v, double theta, double phi)
        {
            return v
                .Rotate(UnitVector3D.XAxis, Angle.FromRadians(theta))
                .Rotate(UnitVector3D.ZAxis, Angle.FromRadians(phi));
        }

        public static Vector3D Rotate(this Vector3D v, Vector3D relativeTo, double theta, double phi)
        {
            return (v - relativeTo).Rotate(theta, phi) + relativeTo;
        }

        public static Vector3D[] RotatePhi(this Vector3D[] vectors, Vector3D relativeTo, double phi)
        {
            var result = new Vector3D[vectors.Length];

            for (int i = 0; i < vectors.Length; i++)
                result[i] = vectors[i].RotatePhi(relativeTo, phi);

            return result;
        }

        public static Vector3D[] Rotate(this Vector3D[] vectors, Vector3D relativeTo, double theta, double phi)
        {
            var result = new Vector3D[vectors.Length];

            for (int i = 0; i < vectors.Length; i++)
                result[i] = vectors[i].Rotate(relativeTo, theta, phi);

            return result;
        }

        public static double Phi(this Vector3D v)
        {
            return Math.Atan2(v.Y, v.X);
        }

        public static double Theta(this Vector3D v)
        {
            return Math.Atan2(v.Z, v.Length);
        }

        public static double Phi(this UnitVector3D v)
        {
            return Math.Atan2(v.Y, v.X);
        }

        public static double Theta(this UnitVector3D v)
        {
            return Math.Atan2(v.Z, v.Length);
        }

        public static string ToBitsString(this uint x, bool doNotPad = true)
        {
            var bitStrings = new string[32];

            int i = 31;

            if (doNotPad)
            {
                while ((i >= 0) && (((x >> i) & 1U) == 0))
                    i--;
            }

            for (; i >= 0; i--)
                bitStrings[31 - i] = ((x >> i) & 1U).ToString();

            return string.Join(string.Empty, bitStrings);
        }

        public static string ToBitsString(this ulong x, bool doNotPad = true)
        {
            var bitStrings = new string[64];

            int i = 63;

            if (doNotPad)
            {
                while ((i >= 0) && (((x >> i) & 1U) == 0))
                    i--;
            }

            for (; i >= 0; i--)
                bitStrings[63 - i] = ((x >> i) & 1U).ToString();

            return string.Join(string.Empty, bitStrings);
        }
    }
}
