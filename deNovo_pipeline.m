path = '~/new_cleaning/data/NDV_sispa_run8/';
load ~/new_cleaning/ndv_r8_files 

adapter = 'GCCGGAGCTCTGCAGATATCG';

samps = length(files):-1:1;

for i=samps
    file = [path, files{i}];
    
     [cloud{i},con{i}, seq{i}]=trinity(file,adapter);
    i
    %[cloud{i},con{i}, seq{i}] = cap3( file, adapter, i );
end
%%
% for i=samps
%     length(con{i})
% end