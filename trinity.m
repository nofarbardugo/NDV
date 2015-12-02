function [ cloud, con, seqs ] = trinity( file, parameters )
% runs trinity to de novo build contig.
% input:
% file = path to input fastq file
% parameters = cleaning parameters
% output:
% cloud = a cloud made of all related trinity sequences combined
% con = consensus sequence of trinity sequnces combined
% seqs = trinity sequences

% remove adapters, trim edges and remove low quality reads
file = clean_fastq(file,'clean',parameters);

% run trinity
unix(sprintf('~/tools/trinityrnaseq-2.1.1/Trinity --seqType fq --single %s --CPU 6 --max_memory 10G',file));

% get trinity results
[~, seqs] = fastaread('trinity_out_dir/Trinity.fasta');
unix('rm -R trinity_out_dir');

% combine trinity results to one sequence
[~, newcastle_con] = fastaread('~/new_cleaning/NDV_full_con.fasta');
[cloud, con] = deNovo_check(seqs, parameters.seq_ln ,...
    parameters.seq_fit_per, newcastle_con);

end

