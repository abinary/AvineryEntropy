function [y] = CatFields(x)

%%
y = [];
fields = fieldnames(x);
for i = 1:numel(fields)
    values = {x.(fields{i})};
    
    all_horizontal = all(cellfun(@(c)(size(c, 1) == 1), values));
    
    if (all_horizontal)
        y.(fields{i}) = horzcat(x.(fields{i}));
    else
        y.(fields{i}) = vertcat(x.(fields{i}));
    end
end

end
