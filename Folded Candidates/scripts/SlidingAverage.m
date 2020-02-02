function [slidingAvg, slidingAvgIndex] = SlidingAverage(x, windowLength, windowStride)

framesIndexes = 1:numel(x, 1);
movingAvgKernel = ones(1, windowLength) * (1/windowLength);

slidingAvg = conv(x, movingAvgKernel, 'valid');
slidingAvg = slidingAvg(1:windowStride:end);

if (nargout >= 2)
    slidingAvgIndex = conv(framesIndexes, movingAvgKernel, 'valid');
    slidingAvgIndex = slidingAvgIndex(1:windowStride:end);
end

end
