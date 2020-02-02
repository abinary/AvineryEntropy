function [y] = VertCatFields(x)

%%
y = [];
fields = fieldnames(x);
for i = 1:numel(fields)
    y.(fields{i}) = vertcat(x.(fields{i}));
end

end
