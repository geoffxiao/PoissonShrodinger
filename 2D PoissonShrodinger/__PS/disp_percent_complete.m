function [flag, display] = disp_percent_complete(curr_loop, tot_loop)

    increment = tot_loop / 10;
    vals = [1 2 3 4 5 6 7 8 9] * increment;
    vals = floor(vals);
    finding = ((curr_loop - vals) == 1);

    if(sum(finding) > 0)
       
        [maxval index] = max(finding);
        flag = true;
        display = index * 10;
        
    else
    
        flag = false;
        display = 0;
    
    end
end

