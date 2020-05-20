using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Avinery.Compression
{
    public class BitArray
    {
        public BitArray(byte[] data = null)
        {
            Data = data;
        }

        public BitArray(int length)
            : this(new byte[length / 8])
        {
        }

        public byte[] Data { get; set; }
        public int Length => Data.Length * 8;

        public BitArray KeepAndSkip(int[] keepSkipPattern)
        {
            int totalToSkip = 0;
            int totalToKeep = 0;

            for (int i = 0; i < keepSkipPattern.Length; i++)
            {
                if (i % 2 == 0)
                    totalToKeep += keepSkipPattern[i];
                else
                    totalToSkip += keepSkipPattern[i];
            }

            double sizeFactor = totalToKeep / (double)(totalToSkip + totalToKeep);

            var newSize = (int)(sizeFactor * Length);

            var result = new BitArray(newSize);

            var currentIndex = 0;
            var currentOutputIndex = 0;
            while (currentIndex < Length)
            {
                for (int i = 0; i < keepSkipPattern.Length; i++)
                {
                    if (i % 2 == 0)
                    {
                        for (int j = 0; j < keepSkipPattern[i]; j++)
                            result[currentOutputIndex + j] = this[currentIndex + j];

                        currentOutputIndex += keepSkipPattern[i];
                    }

                    currentIndex += keepSkipPattern[i]; // skip
                }
            }

            return result;
        }

        public int this[int index]
        {
            get
            {
                var byteIndex = index / 8;
                var bitIndex = index % 8;

                return ((Data[byteIndex] >> bitIndex) & 0x1);
            }

            set
            {
                var byteIndex = index / 8;
                var bitIndex = index % 8;
                var bit = (value != 0 ? 1 : 0) << bitIndex;
                var mask = 0xff ^ (0x1 << bitIndex);

                Data[byteIndex] = (byte)((Data[byteIndex] & mask) | bit);
            }
        }

        public long GetNumNotZero()
        {
            long result = 0;

            for (int i = 0; i < Length; i++)
                if (this[i] != 0)
                    result++;

            return result;
        }

        public long GetNumZero()
        {
            long result = 0;

            for (int i = 0; i < Length; i++)
                if (this[i] == 0)
                    result++;

            return result;
        }

        public BitArray Clone()
        {
            return new BitArray() { Data = (byte[])this.Data.Clone() };
        }
    }
}
