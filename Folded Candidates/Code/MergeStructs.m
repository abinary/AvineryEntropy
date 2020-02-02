function [merged] = MergeStructs(a, b)

%%
merged = a;

fields = fieldnames(b);
for i = 1:numel(fields)
    [merged.(fields{i})] = b.(fields{i});
    %[merged.(fields{i})] = deal({b.(fields{i})});
end

end
