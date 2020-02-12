using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;

namespace Avinery.Sim
{
    public static class Vec2dExtentions
    {
        public static double[] ToDoubleArray(this Vec2d v)
        {
            return new double[2] { v.x, v.y };
        }

        public static double[,] ToDoubleArray(this Vec2d[] v)
        {
            var result = new double[v.Length, 2];

            for (int i = 0; i < v.Length; i++)
            {
                result[i, 0] = v[i].x;
                result[i, 1] = v[i].y;
            }

            return result;
        }

        public static Vec2d ToVec2d(this double[] v)
        {
            var result = new Vec2d();

            result.x = v[0];
            result.y = v[1];

            return result;
        }

        public static Vec2d[] ToVec2dArray(this double[,] v)
        {
            var result = new Vec2d[v.GetLength(0)];

            for (int i = 0; i < result.Length; i++)
                result[i] = new Vec2d() { x = v[i, 0], y = v[i, 1], };

            return result;
        }

        public static Vec2d[] Add(this Vec2d[] vectors, Vec2d v)
        {
            var result = new Vec2d[vectors.Length];

            for (int i = 0; i < vectors.Length; i++)
                result[i] = vectors[i] + v;

            return result;
        }

        public static Vec2d[] Rotate(this Vec2d[] vectors, double angleOfRotation, Vec2d relativeTo)
        {
            var result = new Vec2d[vectors.Length];

            for (int i = 0; i < vectors.Length; i++)
                result[i] = vectors[i].Rotate(angleOfRotation, relativeTo);

            return result;
        }
    }

    [DebuggerDisplay("[{x} , {y}]")]
    public struct Vec2d
    {
        public Vec2d(double x, double y)
        {
            this.x = x;
            this.y = y;
        }

        public static Vec2d FromPolar(double r, double theta)
        {
            return new Vec2d(r * Math.Cos(theta), r * Math.Sin(theta));
        }


        public double x;
        public double y;

        public double Angle { get { return Math.Atan2(y, x); } }
        public double Radius { get { return Norm(); } }
        public Vec2d DirectionVector => Div(Norm());
        public Vec2d OrthogonalVector => new Vec2d() { x = this.y, y = -this.x };

        public static double Dot(Vec2d a, Vec2d b)
        {
            return a.x * b.x + a.y * b.y;
        }

        public double Dot(Vec2d other)
        {
            return this.x * other.x + this.y * other.y;
        }

        public double Norm()
        {
            return Math.Sqrt(this.x * this.x + this.y * this.y);
        }

        public double DistanceSquared(Vec2d from)
        {
            var x = (this.x - from.x);
            var y = (this.y - from.y);
            return (x * x + y * y);
        }

        public double Distance(Vec2d from)
        {
            var x = (this.x - from.x);
            var y = (this.y - from.y);
            return Math.Sqrt(x * x + y * y);
        }

        public static Vec2d operator +(Vec2d a, double b)
        {
            a.x += b;
            a.y += b;
            return a;
        }

        public static Vec2d operator -(Vec2d a, double b)
        {
            a.x -= b;
            a.y -= b;
            return a;
        }

        public static Vec2d operator +(Vec2d a, Vec2d b)
        {
            a.x += b.x;
            a.y += b.y;
            return a;
        }

        public static Vec2d operator -(Vec2d a, Vec2d b)
        {
            a.x -= b.x;
            a.y -= b.y;
            return a;
        }

        public static Vec2d operator *(Vec2d a, Vec2d b)
        {
            a.x *= b.x;
            a.y *= b.y;
            return a;
        }

        public static Vec2d operator /(Vec2d a, Vec2d b)
        {
            a.x /= b.x;
            a.y /= b.y;
            return a;
        }

        public static Vec2d operator %(Vec2d a, Vec2d b)
        {
            a.x %= b.x;
            a.y %= b.y;
            return a;
        }

        public static Vec2d operator *(double b, Vec2d a)
        {
            a.x *= b;
            a.y *= b;
            return a;
        }

        public static Vec2d operator *(Vec2d a, double b)
        {
            a.x *= b;
            a.y *= b;
            return a;
        }

        public static Vec2d operator /(Vec2d a, double b)
        {
            a.x /= b;
            a.y /= b;
            return a;
        }

        public static Vec2d operator %(Vec2d a, double b)
        {
            a.x %= b;
            a.y %= b;
            return a;
        }
        public static Vec2d operator >(Vec2d a, Vec2d b)
        {
            a.x = (a.x > b.x ? 1.0 : 0.0);
            a.y = (a.y > b.y ? 1.0 : 0.0);
            return a;
        }

        public static Vec2d operator >=(Vec2d a, Vec2d b)
        {
            a.x = (a.x >= b.x ? 1.0 : 0.0);
            a.y = (a.y >= b.y ? 1.0 : 0.0);
            return a;
        }

        public static Vec2d operator <(Vec2d a, Vec2d b)
        {
            a.x = (a.x < b.x ? 1.0 : 0.0);
            a.y = (a.y < b.y ? 1.0 : 0.0);
            return a;
        }

        public static Vec2d operator <=(Vec2d a, Vec2d b)
        {
            a.x = (a.x <= b.x ? 1.0 : 0.0);
            a.y = (a.y <= b.y ? 1.0 : 0.0);
            return a;
        }

        public void MultiplyInPlace(double by)
        {
            this.x *= by;
            this.y *= by;
        }

        public void DivInPlace(double by)
        {
            MultiplyInPlace(1.0 / by);
        }

        public Vec2d Multiply(double by)
        {
            var result = new Vec2d();
            result.x = this.x * by;
            result.y = this.y * by;
            return result;
        }

        public Vec2d Div(double by)
        {
            return Multiply(1.0 / by);
        }

        public Vec2d Mod(double by)
        {
            var result = new Vec2d();
            result.x = this.x % by;
            result.y = this.y % by;
            return result;
        }

        public Vec2d Floor()
        {
            return new Vec2d() { x = Math.Floor(x), y = Math.Floor(y)};
        }

        public Vec2d Ceiling()
        {
            return new Vec2d() { x = Math.Ceiling(x), y = Math.Ceiling(y) };
        }

        public Vec2d Round()
        {
            return new Vec2d() { x = Math.Round(x), y = Math.Round(y) };
        }

        public Vec2d Rotate(double angleOfRotation)
        {
            var result = new Vec2d();

            //var angle = Math.Atan2(this.y, this.x);
            //var r = this.Norm();

            //result.x = r*Math.Cos(angle + angleOfRotation);
            //result.y = r*Math.Sin(angle + angleOfRotation);

            var c = Math.Cos(angleOfRotation);
            var s = Math.Sin(angleOfRotation);
            result.x = c*this.x - s*this.y;
            result.y = s*this.x + c*this.y;

            return result;
        }

        public Vec2d Rotate(double angleOfRotation, Vec2d relativeTo)
        {
            var result = (this - relativeTo).Rotate(angleOfRotation) + relativeTo;
            return result;
        }
    }
}
