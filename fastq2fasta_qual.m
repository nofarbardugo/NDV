function [] = fastq2fasta_qual( file_name )
% write fasta and qual files fitting input fastq file.

%faid = fopen(regexprep(file_name,'fastq','fasta'), 'w');
qid = fopen(regexprep(file_name,'fastq','qual'), 'w');
fid = fopen(file_name);

tline= fgets(fid); 
while ischar(tline)
    
    id = tline(2:end-1);
    %getting seq info
    tline= fgets(fid);
    seq = tline(1:end-1);
	tline= fgets(fid);
	tline= fgets(fid);
	quality = tline(1:end-1);
    tline = fgets(fid);

    %fprintf(faid, '>%s\n%s\n', id, seq );
    fprintf(qid, '>%s\n', id );
    for i=1:length(quality)-1
        fprintf(qid, '%d ', quality(i)-33 );
    end
    fprintf(qid, '%d\n', quality(end)-33);
    
end

%fclose(faid);
fclose(qid);
fclose(fid);

end