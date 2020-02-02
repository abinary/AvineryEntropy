function output_txt = WindowDetailsDataCursorFun(src, event_obj)
% ~            Currently not used (empty)
% event_obj    Object containing event data structure
% output_txt   Data cursor text

%pos = get(event_obj,'Position');

eventdata = get(event_obj);
index = get(event_obj, 'DataIndex'); % An undocumented property?

% Get figure
f = event_obj.Target.Parent.Parent;
f.UserData.SelectedIndex = index;

pos = eventdata.Position;

output_txt = sprintf('H: %0.1f S: %0.1f RMSD: %0.1f ', pos);
datacursormode off;
   
end
