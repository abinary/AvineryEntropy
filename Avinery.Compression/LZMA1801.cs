using System;

namespace Avinery.Compression
{
    public static class LZMA1801
    {
        public static int[] SlidingCompressedSize(byte[] data, int stride, int len, int dictionarySizeLog2 = 26)
        {
            if (data.Length < len)
                return new int[0];

            var end = data.Length - len + 1;
            var numWindows = end / stride + 1;
            var results = new int[numWindows];

            try
            {
                var encoder = new lzma1801.LzmaEncoder();
                encoder.DictionarySizeLog2 =
                    dictionarySizeLog2; // Do not go over 26 !!! will consume ~10 times more memory than specified

                for (int idx = 0; idx < numWindows; idx++)
                {
                    var pos = idx * stride;
                    results[idx] = encoder.EncodedLen(data, pos, len);
                }

            }
            catch
            {
            }
            finally
            {
            }

            return results;
        }

        public static long CompressedSize(byte[] data, int dictionarySizeLog2 = 26)
        {
            //var outputBuffer = new byte[(int)(data.Length * 1.2)];
            var encoder = new lzma1801.LzmaEncoder();
            encoder.DictionarySizeLog2 = dictionarySizeLog2; // Do not go over 26 !!! will consume ~10 times more memory than specified

            //var compressedSize = encoder.EncodedLen(data, data.Length, outputBuffer);
            var compressedSize = encoder.EncodedLen(data, 0, data.Length);

            return compressedSize;
        }
    }
}