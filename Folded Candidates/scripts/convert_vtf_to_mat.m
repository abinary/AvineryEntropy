function convert_vtf_to_mat(filename)

f = fopen(filename, 'r');

l = fgetl(f);
while (~strcmp(l, 'timestep ordered'))
    l = fgetl(f);
end

frames = [];
i = 0;
while (~feof(f))
    i = i + 1;
    numbers = textscan(f, '%f');
    numbers = numbers{1};
    
    numbers = reshape(numbers, [3, numel(numbers)/3])';
    
    if (i == 1)
        frames = reshape(numbers, [1 size(numbers)]);
    end
    
    % Cause preallocation of more frames
    if (mod(i, 5000) == 1)
        frames(ceil(i / 5000) * 5000, end, end) = 0;
    end
    
    frames(i, :, :) = numbers;
    
    l = fgetl(f);
end

fclose(f);

frames  = frames(1:i, :, :);
save([filename '.mat'], 'frames');


end
