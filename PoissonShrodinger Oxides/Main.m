%%
clear all;
clc

% 10 u.c. LAO on STO at different Temperatures
start_time = clock;
save('start_time');

%%

clear all; clc;
Potential_Guess_LAO_STOSub_300
Start_PoissonShrodinger

%%
clear all;
finish_time = clock;
save('finish_time1');
%%

clear all; clc;
Potential_Guess_LAO_STOSub_4
Start_PoissonShrodinger

%%
clear all;
finish_time = clock;
save('finish_time2');
