function [final_file] = clean_fastq( file_name, out_name, parameters )

% trim_file = regexprep(file_name,'\.fastq','_qualTrm.fastq');
% adapters_file = regexprep(file_name,'\.fastq','_noAd.fastq');
final_file = regexprep(file_name,'\.fastq',sprintf('_%s.fastq',out_name));

% unix( sprintf('~/tools/cutadapt-1.8.3/bin/cutadapt -a %s -o %s %s', adapter, adapters_file, file_name))
% unix( sprintf('fastx_clipper -a %s -i %s -o %s',adapter, file_name, adapters_file));
% unix(sprintf('fastx_trimmer -f 12 -i %s -o %s',adapters_file,final_file));
% trim_and_filter( adapters_file, final_file, 25, 0,10, 20, 20 )
% unix(sprintf('fastq_quality_filter -q 30 -i %s -o %s',adapters_file,final_file));

% unix( sprintf('~/tools/cutadapt-1.8.3/bin/cutadapt -a %s -o %s %s', adapter, adapters_file, file_name))
% trim_and_filter( adapters_file, final_file, 25, 0,10, 20, 20 )
% unix(sprintf('rm %s',adapters_file));

trim_and_filter( file_name, final_file, parameters.len, ...
    parameters.trim_start,parameters.trim_end, ...
    parameters.minPer, parameters.minVal )


end

