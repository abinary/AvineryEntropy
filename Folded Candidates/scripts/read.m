
f = fopen('nug2_0.vtf', 'r');

l = fgetl(f);
while (~strcmp(l, 'timestep ordered'))
    l = fgetl(f);
end

frames = [];

i = 0;
while (~feof(f))
    i = i + 1
    numbers = textscan(f, '%f');
    numbers = numbers{1};
    
    numbers = reshape(numbers, [3, numel(numbers)/3])';
    
    frames(end + 1, :, :) = numbers;
    
    l = fgetl(f);
end

fclose(f);

save('nug2_0.vtf.mat', 'frames');
