using MathNet.Numerics.Random;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Avinery.Sim
{
    public class MCBase
    {
        public MCBase()
        {
            var r = new Random();
            RandomSeed = r.Next();
        }

        private int _randomSeed;

#if true
        protected MersenneTwister _random;
        public int RandomSeed
        {
            get { return _randomSeed; }
            set { _randomSeed = value; _random = new MersenneTwister(_randomSeed); }
        }

#else
        protected Random _random;
        public int RandomSeed
        {
            get { return _randomSeed; }
            set { _randomSeed = value; _random = new Random(_randomSeed); }
        }
#endif
    }
}
