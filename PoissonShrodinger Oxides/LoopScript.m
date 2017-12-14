clear all;
clc;

dEc12s = [0.1 : 0.2 : 1.5]; % 2 is above 3 by this much

for loop_counter_variable = 1 : numel(dEc12s)
    
    dEc12 = dEc12s(loop_counter_variable);
    FILE_NAME = strcat('CSO_BLSO_Offset_',num2str(dEc12s(loop_counter_variable)));
    Start_PoissonShrodinger
    clearvars -except loop_counter_variable dEc12s dEc12 FILE_NAME
     
end