// lzma1801.h

#pragma once

//using namespace System;
using namespace System::Runtime::InteropServices;

#include "C/LzmaLib.h"
#include "C/LzmaEnc.h"
#include "C/Alloc.h"

typedef unsigned char byte;
typedef byte byte_array[];
typedef byte* byte_ptr;

namespace lzma1801 {

	public ref class LzmaInternalStatistics
	{
	public:
		int numOfMatches;

		long sumOfMatchLength;
		long sumOfMatchLengthSquare;

		long sumOfMatchDistance;
		long sumOfMatchDistanceSquare;

		double sumOfMatchLogDistance;
		double sumOfMatchLogDistanceSquare;

		long maxMatchLength;
		long minMatchLength;
		long maxMatchDistance;
		long minMatchDistance;
	};

	public ref class LzmaEncoder
	{
	public:

		int Level = 9;
		int DictionarySizeLog2 = 26;
		int LiteralContextBits = 3; // 0-8
		int LiteralPosBits = 0; // 0-4
		int NumOfPosBits = 2; // 0-4
		int NumOfFastBytes = 64;
		int BTMode = 1;
		int NumOfHashBytes = 4;

		long Encode(array<System::Byte>^ input, int inputLength, array<System::Byte>^ output, int maxOutputLength)
		{
			pin_ptr<System::Byte> inputPtr = &input[0];
			pin_ptr<System::Byte> outputPtr = &output[0];
			return Encode(inputPtr, inputLength, outputPtr, maxOutputLength);
		}

		long Encode(array<System::Byte>^ input, int startPos, int inputLength, array<System::Byte>^ output, int outputPos, int maxOutputLength)
		{
			pin_ptr<System::Byte> inputPtr = &input[startPos];
			pin_ptr<System::Byte> outputPtr = &output[outputPos];
			return Encode(inputPtr, inputLength, outputPtr, maxOutputLength);
		}

		long Encode(array<System::Byte>^ input, int inputLength, array<System::Byte>^ output)
		{
			pin_ptr<System::Byte> inputPtr = &input[0];
			pin_ptr<System::Byte> outputPtr = &output[0];
			return Encode(inputPtr, inputLength, outputPtr, output->Length);
		}

		long Encode(array<System::Byte>^ input, array<System::Byte>^ output)
		{
			pin_ptr<System::Byte> inputPtr = &input[0];
			pin_ptr<System::Byte> outputPtr = &output[0];
			return Encode(inputPtr, input->Length, outputPtr, output->Length);
		}

		long Encode(array<System::Byte>^ input)
		{
			array<System::Byte>^ output = gcnew array<System::Byte>(input->Length + input->Length / 2);
			return Encode(input, output);
		}

		long _Encode(byte_ptr input, int inputLength, byte_ptr output, int maxOutputLength)
		{
			size_t outputLength = maxOutputLength;
			size_t outPropsSize = 5;
			byte outProps[5];

			auto result = LzmaCompress(output, &outputLength, input, inputLength, outProps, &outPropsSize,
				Level, 1 << DictionarySizeLog2, LiteralContextBits, LiteralPosBits, NumOfPosBits, NumOfFastBytes, 1);

			return outputLength;
		}

		long Encode(byte_ptr input, int inputLength, byte_ptr output, int maxOutputLength)
		{
			size_t outputLength = maxOutputLength;
			size_t outPropsSize = 5;
			byte outProps[5];

			CLzmaEncProps props;
			LzmaEncProps_Init(&props);
			props.level = Level;
			props.dictSize = 1 << DictionarySizeLog2;
			props.lc = LiteralContextBits;
			props.lp = LiteralPosBits;
			props.pb = NumOfPosBits;
			props.fb = NumOfFastBytes;
			props.btMode = BTMode;
			props.numHashBytes = NumOfHashBytes;
			props.numThreads = 1;

			ICompressProgress progress;

			auto progressDelegate = gcnew ProgressInternalDelegate(this, &LzmaEncoder::Progress);

			// When the code is optimized in "release" build, the delegate seems like it is unused after the
			// conversion to a pointer and is sometimes collected by the garbage collector. To avoid that,
			// the GCHandle is used.
			auto progressDelegateHandle = GCHandle::Alloc(progressDelegate);

			auto progressDelegatePtr = Marshal::GetFunctionPointerForDelegate(progressDelegate);
			progress.Progress = (SRes(*)(const ICompressProgress * p, UInt64 inSize, UInt64 outSize))(void *) progressDelegatePtr;

			//SetMaxUnpackedBlockSize(65536U);

			//auto resultCode = LzmaEncode(output, &outputLength, input, inputLength, &props, outProps, &outPropsSize, 0,
			//	NULL, &g_Alloc, &g_Alloc);
			auto resultCode = LzmaEncode(output, &outputLength, input, inputLength, &props, outProps, &outPropsSize, 0,
				&progress, &g_Alloc, &g_Alloc);

			progressDelegateHandle.Free();

			return outputLength;
		}

		long EncodedLen(array<System::Byte>^ inputArray, int inputPos, int inputLength)
		{
			pin_ptr<System::Byte> inputPtr = &inputArray[0];
			return EncodedLen(inputPtr, inputPos, inputLength);
		}

		//long EncodedLen(GCHandle input_fixed, int inputPos, int inputLength)
		//{
		//	pin_ptr<System::Byte> inputPtr = (byte_ptr)GCHandle::ToIntPtr(input_fixed).ToInt64();
		//	//byte_ptr input = (byte_ptr)GCHandle::ToIntPtr(input_fixed).ToInt64();
		//	return EncodedLen(inputPtr, inputPos, inputLength);
		//}

		void EncodedLen(array<System::Byte>^ inputArray, int inputPos, int inputLength, array<System::Int64>^ outputArray, int outputPos)
		{
			pin_ptr<System::Byte> inputPtr = &inputArray[0];
			outputArray[outputPos] = EncodedLen(inputPtr, inputPos, inputLength);
		}

		long EncodedLen(byte_ptr input, int inputPos, int inputLength)
		{
			size_t outputLength = 1<<31;
			size_t outPropsSize = 5;
			byte outProps[5];

			CLzmaEncProps props;
			LzmaEncProps_Init(&props);
			props.level = Level;
			props.dictSize = 1 << DictionarySizeLog2;
			props.lc = LiteralContextBits;
			props.lp = LiteralPosBits;
			props.pb = NumOfPosBits;
			props.fb = NumOfFastBytes;
			props.btMode = BTMode;
			props.numHashBytes = NumOfHashBytes;
			props.numThreads = 1;

			ICompressProgress progress;

			auto progressDelegate = gcnew ProgressInternalDelegate(this, &LzmaEncoder::Progress);

			// When the code is optimized in "release" build, the delegate seems like it is unused after the
			// conversion to a pointer and is sometimes collected by the garbage collector. To avoid that,
			// the GCHandle is used.
			auto progressDelegateHandle = GCHandle::Alloc(progressDelegate);

			auto progressDelegatePtr = Marshal::GetFunctionPointerForDelegate(progressDelegate);
			progress.Progress = (SRes(*)(const ICompressProgress * p, UInt64 inSize, UInt64 outSize))(void*) progressDelegatePtr;

			//SetMaxUnpackedBlockSize(65536U);

			//auto resultCode = LzmaEncode(output, &outputLength, input, inputLength, &props, outProps, &outPropsSize, 0,
			//	NULL, &g_Alloc, &g_Alloc);
			auto resultCode = LzmaEncode(NULL, &outputLength, input + inputPos, inputLength, &props, outProps, &outPropsSize, 0,
				&progress, &g_Alloc, &g_Alloc);

			progressDelegateHandle.Free();

			return outputLength;
		}

		LzmaInternalStatistics^ LastStatistics;

		delegate void ProgressDelegate(UInt64 inSize, UInt64 outSize);
		event ProgressDelegate^ OnProgress;

	private:
		delegate SRes ProgressInternalDelegate(void *p, UInt64 inSize, UInt64 outSize);
		SRes Progress(void *p, UInt64 inSize, UInt64 outSize)
		{
			try
			{
				OnProgress(inSize, outSize);
			}
			catch (System::Exception^ ex)
			{
			}

			return SZ_OK;
		}
	};
}
