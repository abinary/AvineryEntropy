﻿
/*
	DO NOT EDIT THIS FILE!!! THIS FILE IS AUTO-GENERATED...
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;

namespace Avinery.Sim
{
    [Serializable]
    public struct Tuple2<T1, T2>
    {
        public T1 Item1;
        public T2 Item2;

        public Tuple2(T1 item1, T2 item2)
        {
            this.Item1 = item1;
            this.Item2 = item2;
        }

        public override bool Equals(object obj)
        {
            if (!(obj is Tuple2<T1, T2>))
                return false;

            var otherTuple = (Tuple2<T1, T2>)obj;

            return
                this.Item1.Equals(otherTuple.Item1) &&
                this.Item2.Equals(otherTuple.Item2);
        }

        public override int GetHashCode()
        {
            return
                this.Item1.GetHashCode() ^
                this.Item2.GetHashCode();
        }
    }

    public struct Tuple3<T1, T2, T3>
    {
        public T1 Item1;
        public T2 Item2;
        public T3 Item3;

        public Tuple3(T1 item1, T2 item2, T3 item3)
        {
            this.Item1 = item1;
            this.Item2 = item2;
            this.Item3 = item3;
        }

        public override bool Equals(object obj)
        {
            if (!(obj is Tuple3<T1, T2, T3>))
                return false;

            var otherTuple = (Tuple3<T1, T2, T3>)obj;

            return
                this.Item1.Equals(otherTuple.Item1) &&
                this.Item2.Equals(otherTuple.Item2) &&
                this.Item3.Equals(otherTuple.Item3);
        }

        public override int GetHashCode()
        {
            return
                this.Item1.GetHashCode() ^
                this.Item2.GetHashCode() ^
                this.Item3.GetHashCode();
        }
    }

    public struct Tuple4<T1, T2, T3, T4>
    {
        public T1 Item1;
        public T2 Item2;
        public T3 Item3;
        public T4 Item4;

        public Tuple4(T1 item1, T2 item2, T3 item3, T4 item4)
        {
            this.Item1 = item1;
            this.Item2 = item2;
            this.Item3 = item3;
            this.Item4 = item4;
        }

        public override bool Equals(object obj)
        {
            if (!(obj is Tuple4<T1, T2, T3, T4>))
                return false;

            var otherTuple = (Tuple4<T1, T2, T3, T4>)obj;

            return
                this.Item1.Equals(otherTuple.Item1) &&
                this.Item2.Equals(otherTuple.Item2) &&
                this.Item3.Equals(otherTuple.Item3) &&
                this.Item4.Equals(otherTuple.Item4);
        }

        public override int GetHashCode()
        {
            return
                this.Item1.GetHashCode() ^
                this.Item2.GetHashCode() ^
                this.Item3.GetHashCode() ^
                this.Item4.GetHashCode();
        }
    }

    public struct Tuple5<T1, T2, T3, T4, T5>
    {
        public T1 Item1;
        public T2 Item2;
        public T3 Item3;
        public T4 Item4;
        public T5 Item5;

        public Tuple5(T1 item1, T2 item2, T3 item3, T4 item4, T5 item5)
        {
            this.Item1 = item1;
            this.Item2 = item2;
            this.Item3 = item3;
            this.Item4 = item4;
            this.Item5 = item5;
        }

        public override bool Equals(object obj)
        {
            if (!(obj is Tuple5<T1, T2, T3, T4, T5>))
                return false;

            var otherTuple = (Tuple5<T1, T2, T3, T4, T5>)obj;

            return
                this.Item1.Equals(otherTuple.Item1) &&
                this.Item2.Equals(otherTuple.Item2) &&
                this.Item3.Equals(otherTuple.Item3) &&
                this.Item4.Equals(otherTuple.Item4) &&
                this.Item5.Equals(otherTuple.Item5);
        }

        public override int GetHashCode()
        {
            return
                this.Item1.GetHashCode() ^
                this.Item2.GetHashCode() ^
                this.Item3.GetHashCode() ^
                this.Item4.GetHashCode() ^
                this.Item5.GetHashCode();
        }
    }
}