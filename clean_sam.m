
function [ pos_read_count, insert_count ] = clean_sam( file_name, o_file_name, source_seq_len )
% takes sam file and writes good aligned reads to fastq file.

pos_read_count = zeros(1,source_seq_len);
insert_count = zeros(1,source_seq_len);

fid = fopen(file_name);
oid = fopen(o_file_name, 'w');

tline= fgets(fid);
while ischar(tline)
     %not header of sam file
    if (tline(1) ~= '@')

        %split line into aligment info 
        alInfo = regexp(tline,'\t','split');
        id = alInfo{1};
        pos = sscanf(alInfo{4},'%d');
        mapQuality = sscanf(alInfo{5},'%d');
        cigar = alInfo{6};
        seq = alInfo{10};
        quality = alInfo{11};
        
        % if alignment is significant
        if mapQuality > 0 
                seq = cigar2align({seq},{cigar});
                % if there is an insert in the read
                if sum(cigar=='I') ~= 0
                    insert_count(pos) = insert_count(pos)+1;
                end
                pos_read_count(pos) = pos_read_count(pos)+1;
                fprintf(oid, '@%.4d_%s\n%s\n+\n%s\n', pos, id, seq, quality );
        end
    end
    tline = fgets(fid);
end

fclose(fid);
fclose(oid);

end