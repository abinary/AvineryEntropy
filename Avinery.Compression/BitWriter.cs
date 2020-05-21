using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization.Formatters.Binary;
using System.Text;

namespace Avinery.Compression
{
    public class BitWriter : BinaryWriter
    {
        public BitWriter(Stream stream)
            : base(stream)
        {
        }

        private ulong _data = 0;
        private int _bitsPartialFillLength = 0;

        int _readBuffer = 0;
        int _bitsInReadBuffer = 0;

        public long Length { get { return BaseStream.Length * 8; } }

        public int BitsWaiting { get { return _bitsPartialFillLength; } }

        public void WriteBits(ulong data, int numOfBits)
        {
            if (_bitsPartialFillLength + numOfBits > 64)
            {
                var numOfWrittenBits = 64 - _bitsPartialFillLength;
                var numOfBitsLeftToWrite = numOfBits - numOfWrittenBits;
                WriteBits(data, numOfWrittenBits);
                WriteBits(data >> numOfWrittenBits, numOfBitsLeftToWrite);
                return;
            }

            var mask = ulong.MaxValue >> (64 - numOfBits);
            _data |= ((data & mask) << _bitsPartialFillLength);

            _bitsPartialFillLength += numOfBits;

            if (_bitsPartialFillLength == 64)
            {
                Write(_data);
                _data = 0;
                _bitsPartialFillLength = 0;
            }
        }

        public void WriteBits(ulong[] data, int numOfBitsEach)
        {
            for (int i = 0; i < data.Length; i++)
                WriteBits(data[i], numOfBitsEach);
        }

        public void FlushBits()
        {
            Write(_data);
            _bitsPartialFillLength = 0;
        }

        bool _wasClosed = false;
        public override void Close()
        {
            if (_bitsPartialFillLength > 0)
                Write(_data);

            _data = 0;
            _bitsPartialFillLength = 0;

            _wasClosed = true;
            base.Close();
        }

        protected override void Dispose(bool disposing)
        {
            if (!_wasClosed)
                Close();
            base.Dispose(disposing);
        }
    }

    public class BitReader : BinaryReader
    {
        public BitReader(Stream stream)
            : base(stream)
        {
        }

        int _readBuffer = 0;
        int _bitsInReadBuffer = 0;

        public long Length { get { return BaseStream.Length * 8; } }

        public ulong ReadBits(int numOfBits)
        {
            ulong result = 0;
            int bitsInResult = 0;

            if (_bitsInReadBuffer == 0)
            {
                _readBuffer = BaseStream.ReadByte();
                _bitsInReadBuffer = 8;
            }

            while (numOfBits > _bitsInReadBuffer)
            {
                var mask = ~(uint.MaxValue << _bitsInReadBuffer);

                result |= (mask & (uint)_readBuffer) << bitsInResult;
                bitsInResult += _bitsInReadBuffer;

                numOfBits -= _bitsInReadBuffer;
                _bitsInReadBuffer = 0;

                // TODO: Some indication if more bits have been requested than there are in the stream
                _readBuffer = BaseStream.ReadByte();
                _bitsInReadBuffer = 8;
            }

            //if (numOfBits <= _bitsInReadBuffer)
            {
                var mask = ~(uint.MaxValue << numOfBits);

                result |= ((mask & (uint)_readBuffer) << bitsInResult);

                _readBuffer >>= numOfBits;

                //numOfBits = 0;
                _bitsInReadBuffer -= numOfBits;
            }

            return result;
        }

        public ulong[] ReadBits(int numOfBits, int numberOfTimes)
        {
            var result = new ulong[numberOfTimes];

            for (int i = 0; i < result.Length; i++)
                result[i] = ReadBits(numOfBits);

            return result;
        }
    }
}
