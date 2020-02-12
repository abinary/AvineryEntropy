

NET.addAssembly([cd '\bin\x64\Release\netstandard2.0\Avinery.Compression.dll']);

data = uint8(randi(256, [1 1e6]));

x = double(Avinery.Compression.LZMA.CompressedSize(data))

stride = 512;
window_len = 1024;
y = double(Avinery.Compression.LZMA.SlidingCompressedSize(data, stride, window_len))
