function [ con, cloud, insert_count ] = read2seq( in_file, db, source_seq, parms )
% gets fastq file of reads and bowtie database. clean the reads
% alignes using bowtie and makes cloud and consensus. 
% in_file = path to input fastq file
% db = path to bowtie database
% sourse_seq = reference sequence for cloud size
% parms = cleaning parameters.
% con = consensus sequnce
% cloud = nucleotide frequncies cloud
% insert_count = number if insert in positions over the genome, for basic
% check only!

% making names for new files
sam_file = regexprep(in_file,'\.fastq','.sam');
final_clean = regexprep(in_file,'\.fastq','.cleanf_fastq');

% cleanign
clean_file = clean_fastq(in_file,'clean',parms);

% aligning
% unix(sprintf('/export/peptibase/kleintz/bowtie-1.0.0/bowtie -e 500 -q -n 3 --maxbts 1300 --tryhard --best -p1 %s %s %s --sam', db, initial_clean, sam_file));
% unix(sprintf('~/tools/bowtie2-2.2.4/bowtie2 --score-min L,-0.4,-0.4 -x %s %s %s --un-conc %s', db, initial_clean, sam_file, [initial_clean '_noAl.fasta']));
% unix(sprintf('~/tools/bowtie2-2.2.4/bowtie2 --local --dpad 60 --gbar 1 --rdg 3,1 -x %s %s %s --un-conc %s', db, initial_clean, sam_file, [initial_clean '_noAl.fasta']));
% unix(sprintf('~/tools/bowtie2-2.2.4/bowtie2 -x %s %s %s --un-conc %s', db, initial_clean, sam_file, [initial_clean '_noAl.fasta']));
unix(sprintf(parms.al_cmd, db, clean_file, sam_file));

source_seq_len = length(source_seq);

% high aligment quality reads of sam to fastq file
[ ~, insert_count ] = clean_sam( sam_file, final_clean, source_seq_len );

%making cloud
[ cloud ] = cleanfastq_to_cloud( final_clean, source_seq_len );

con = getConsensus(cloud);

coverage_fig( cloud, con, in_file )

%removing remeining files
% unix(sprintf('rm -f %s %s %s',clean_file, sam_file, final_clean ));
%unix(sprintf('rm -f %s %s',clean_file, sam_file ));
end
