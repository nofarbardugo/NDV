% do novo building cotings from files 

path = 'data/NDV_sispa_run8/';
load ndv_r8_files 


% % adapter sequence to remove
% parameters.adapter = 'GCCGGAGCTCTGCAGATATCG';
% minimal read length
parameters.len = 30;
% num of nucletides in read start to remove
parameters.trim_start = 20;
% num of nucletides in read end to remove
parameters.trim_end = 10;
% minimum % of low quality nucletides to allow in read
parameters.minPer = 20;
% minimum phred score to consider high-qiality.
parameters.minVal = 53;
% minimal sequence length to use in deNovo_check
parameters.seq_ln = 300;
% minimal % of fit to known sequnce to use in deNovo_check
parameters.seq_fit_per = 0.15;


samps = 1:length(files);

for i=samps
    file = [path, files{i}];
    [cloud{i},con{i}, seq{i}]=trinity(file,parameters);
end

%%
for i=samps
    [~, full_con{i}] = tcloud2con(tcloud{i});
end
save r8_trinity_cons full_con


