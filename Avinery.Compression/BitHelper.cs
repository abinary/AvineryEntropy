using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Avinery.Compression
{
    public static class BitHelper
    {
        public static byte[] CondenseBits(byte[] data, int numOfBitsEach)
        {
            var output_length = (data.Length * numOfBitsEach + 8 - 1) / 8;
            using (var ms = new MemoryStream(output_length))
            {
                var bw = new BitWriter(ms);
                for (int i = 0; i < data.Length; i++)
                    bw.WriteBits(data[i], numOfBitsEach);

                return ms.ToArray();
            }
        }

        public static byte[] ExpandBits(byte[] data)
        {
            var ba = new BitArray(data);
            var output_length = data.Length * 8;
            var output = new byte[output_length];

            for (int i = 0; i < data.Length; i++)
            {
                output[i * 8 + 0] = (byte)(data[i] & (1U << 0));
                output[i * 8 + 1] = (byte)(data[i] & (1U << 1));
                output[i * 8 + 2] = (byte)(data[i] & (1U << 2));
                output[i * 8 + 3] = (byte)(data[i] & (1U << 3));
                output[i * 8 + 4] = (byte)(data[i] & (1U << 4));
                output[i * 8 + 5] = (byte)(data[i] & (1U << 5));
                output[i * 8 + 6] = (byte)(data[i] & (1U << 6));
                output[i * 8 + 7] = (byte)(data[i] & (1U << 7));
            }

            return output;
        }
    }
}
