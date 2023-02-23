%% 
table = ['a','b','c','f','k','z'];
convert = ["he",'lo','alb','vz','lf','is'];

ex = ['ck'];

function chart = replavechr(table, convert)
    while true
        result = [];
        for i = 1:length(ex)
            temp = ex(i);
            idx = find(table == temp);
            result = [result, convert(idx)];
        end
        % 
        chart = [];
        for ii = 1:length(result)
            chart = [chart, result{ii}];
        end
        % check
        for j = 1:length(chart)
            
            
            %
            if isempty(find(table == string(chart)) % 바꿀 char no
                answer = chart;
            else
                ex = chart;
            end
        end               
    end
    
end

       



