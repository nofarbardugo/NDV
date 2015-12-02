function [ cloud, con, seq ] = cap3( file, adapter,num )
% runs cap3 to de novo build contig.
% input:
% file = path to input fastq file
% adapter = adapter sequence to remove
% output:
% cloud = a cloud made of all related sequences combined
% con = consensus sequence of sequnces combined
% seqs = cap3 sequences

% remove adapters, trim edges and remove low quality reads
file = clean_fastq(file,'clean',adapter);

% run cap3 
a_file = file;
a_file(end)='a';

unix(sprintf('fastq_to_fasta -i %s -o %s',file,a_file));
fastq2fasta_qual( file );

unix(sprintf('~/tools/CAP3/cap3 %s > cap3_res/res%d.txt',a_file,num));

% combine results to one sequence
[~, seq] = fastaread([a_file '.cap.contigs']);
[~, newcastle_con] = fastaread('~/new_cleaning/NDV_full_con.fasta');
[cloud, con] = deNovo_check(seq, 100, 0.15, newcastle_con);


end

