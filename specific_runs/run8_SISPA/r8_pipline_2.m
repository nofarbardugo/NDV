
search_path = 'data/NDV_sispa_run8/*.fastq';
path = 'data/NDV_sispa_run8/';

% setting cleanig parameters

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
% alignment command for bowtie
%al_cmd = '~/tools/bowtie2-2.2.4/bowtie2 --local --dpad 40 --gbar 1 --rdg 5,1 -x %s %s %s';
parameters.al_cmd = '~/tools/bowtie2-2.2.4/bowtie2 --dpad 40 --gbar 1 --rdg 5,1 -x %s %s %s';

% % getting files
% temp=dir(search_path);
% files={temp.name}';
%%

% making newcastle and lasota bowtie databases

db_file{1} = 'db_files/VH_bowtie_index/VH.fasta';
db{1} = 'VH_bowtie2_index/laSota_israel';

db_file{2} = '~/data/NDV/NDV_laSota_28_10_14_noFar.al.fasta';
db{2} = 'VH_bowtie2_index/newcastle_all_v2_noFar';

[~, newcastle_con] = fastaread('db_files/NDV_full_con.fasta');
load('specific_runs/run8_SISPA/r8_trinity_cons.mat')

temp_seq{1} = newcastle_con;
temp_id{1} = 'newcastle consensus';
temp_id{2} = 'de novo sequence';

db_file{3} = 'db_files/temp.fasta';
db{3} = 'VH_bowtie2_index/temp';

used_db_index = 3;

% building bowtie index
% unix(sprintf('~/tools/bowtie2-2.2.4/bowtie2-build %s %s',db_file{used_db_index}, db{used_db_index}));

%%
[~, newcastle_con] = fastaread('db_files/NDV_full_con.fasta');
[~,laSota_con] = fastaread('db_files/VH.fasta');

load ndv_r8_files
load r8_trinity_cons

samps = 1:length(files);

%runing for each file
for i=samps;
    % setting bowtie database to use on file. only tha last file 
    % is lasota, so the dataset will be different
    used_db = db{used_db_index};
    source_seq = newcastle_con;
    
    % making db to fit consensus obtained by de novo
    unix(sprintf('rm %s',db_file{used_db_index}));
    temp_seq{2} = full_con{i};
    fastawrite(db_file{used_db_index},temp_id,temp_seq);
    unix(sprintf('~/tools/bowtie2-2.2.4/bowtie2-build %s %s',db_file{used_db_index}, db{used_db_index}));

%     if i==2
%         used_db = db{2};
%         source_seq = laSota_con;
%     end
    file = [path, files{i}];
    [ consensuses{i}, clouds{i} ] = read2seq(file, used_db, source_seq, parameters );
    
end

%%
for i=samps
    % clouds are larger than refernce, in case of gap. this is the fix
    while consensuses{i}(end) =='N' &&  length(consensuses{i}) > 15192
        consensuses{i}(end) = [];
    end
end

%% ploting
for i=samps
    ids{i} = sprintf('NDV %d',i);
    coverage_fig( clouds{i}, consensuses{i}, ids{i})
end

% %%
% for i=samps
%     noReads{i} = find(consensuses{i}=='N');
%     gaps{i} = find(consensuses{i}=='-');
% end

%% removing gaps and checking coverage
for i=samps
     consensuses{i}(consensuses{i}=='-')=[];
     covs(i) = sum(consensuses{i}~='N')/length(consensuses{i});
end

%% writing result to file
load r8_names
% fastawrite('res/NDV_run8_sequences.fasta',ids,consensuses);
% unix( '~/tools/clustalo -i NDV_run6_sequences.fasta -o NDV_run6_sequences.al.fasta -v');
