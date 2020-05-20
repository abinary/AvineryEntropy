using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Avinery.Compression
{
    public static class Layout
    {
        public static byte[] ApplyByteMask(byte[] input, byte mask)
        {
            var result = new byte[input.Length];
            for (int i = 0; i < input.Length; i++)
                result[i] = (byte)(input[i] & mask);

            return result;
        }

        public static int[] Range(int count)
        {
            var result = new int[count];
            for (int i = 0; i < count; i++)
                result[i] = i;

            return result;
        }

        public static BitArray ReorderBits(BitArray original, int[] order)
        {
            var resultArray = new BitArray(original.Length);

            var chunkSizeInBits = order.Length;

            var numOfChunks = original.Length / chunkSizeInBits;
            for (int chunkIndex = 0; chunkIndex < numOfChunks; chunkIndex++)
            {
                var chunkOffset = chunkIndex * chunkSizeInBits;

                for (int i = 0; i < chunkSizeInBits; i++)
                    resultArray[chunkOffset + i] = original[chunkOffset + order[i]];
            }

            return resultArray;
        }

        public static byte[] ReorderBytes(byte[] original, int[] order)
        {
            var resultArray = new byte[original.Length];
            var chunkLength = order.Length;

            var numOfChunks = original.Length / chunkLength;
            for (int chunkIndex = 0; chunkIndex < numOfChunks; chunkIndex++)
            {
                var chunkOffset = chunkIndex * chunkLength;

                for (int i = 0; i < chunkLength; i++)
                    resultArray[chunkOffset + i] = original[chunkOffset + order[i]];
            }

            return resultArray;
        }
    }
}
