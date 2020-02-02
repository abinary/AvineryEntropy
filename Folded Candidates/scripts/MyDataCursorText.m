function output_txt = MyDataCursorText(~, event_obj)
% ~            Currently not used (empty)
% event_obj    Object containing event data structure
% output_txt   Data cursor text
eventdata = get(event_obj);
pos = eventdata.Position;
output_txt = num2str(pos);
