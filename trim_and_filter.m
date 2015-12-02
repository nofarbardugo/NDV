function [] = trim_and_filter( in_file, out_file, len, trim_start,trim_end, minPer, minVal )
% get fastq file. trims, remove sequnces that are too short
% or have too low quality.
% write result to out_file.
% parameters =  cleaning parameters, includes:
% len = minimal read length
% trim_start = num of nucletides in read start to remove
% trim_end = num of nucletides in read end to remove
% minPer = minimum % of low quality nucletides to allow in read
% minVal = minimum phred score to consider high-qiality.

oid = fopen(out_file, 'w');
fid = fopen(in_file);

tline= fgets(fid); 
while ischar(tline)
    
	id = tline(1:end-1);
    %getting seq info
    tline = fgets(fid);
    seq = tline(trim_start+1:end-(trim_end+1));
	tline = fgets(fid);
	tline = fgets(fid);
	quality = tline(trim_start+1:end-(trim_end+1));
    tline = fgets(fid);

    % if read is good
    if length( seq ) >= len && checkSeq( quality, minPer, minVal ) 
 
        fprintf(oid, '%s\n%s\n+\n%s\n', id, seq, quality);
        
    end
	
end

fclose(oid);
fclose(fid);


end

function [ res ] = checkSeq( quality, minPer, minVal )
% checks if there are more low qulity(lower than minVal)
% reads in seq than minPer%.
% 1 if there aren't.(=1 if clean)

len = length(quality);
quality( uint8(quality) < minVal ) = [];

res = ( len - length(quality) < (minPer/100)*len );

end

