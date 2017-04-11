function colnum=find_column_number(hdrs,field)

colnum = 0;
for i=1:length(hdrs)
   if strcmpi(hdrs{i},field),
       colnum = i;
   end
end