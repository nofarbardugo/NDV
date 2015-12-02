
function [ cloud1, cloud2, cloud3 ] = cleanfastq_to_cloud( file_name, template_len )
%gets aligned clean fastq file and creates 
%clouds of necleotides, nucleotides pairs and codons.

nucs = {'A', 'C', 'G', 'T', 'N', '-'};
nucPairs = {'AA', 'AC', 'AG', 'AT', 'AN', 'A-',  ...
		  'CA', 'CC', 'CG', 'CT', 'CN', 'C-',  ...
		  'GA', 'GC', 'GG', 'GT', 'GN', 'G-',  ...
		  'TA', 'TC', 'TG', 'TT', 'TN', 'T-',  ...
		  'NA', 'NC', 'NG', 'NT', 'NN', 'N-',  ...
		  '-A', '-C', '-G', '-T', '-N', '--'};
codons = {'AAA', 'AAC', 'AAG', 'AAT', 'AAN', 'AA-', 'ACA', 'ACC', 'ACG', 'ACT', 'ACN', 'AC-', 'AGA', 'AGC', 'AGG', 'AGT', 'AGN', 'AG-', 'ATA', 'ATC', 'ATG', 'ATT', 'ATN', 'AT-', 'ANA', 'ANC', 'ANG', 'ANT', 'ANN', 'AN-', 'A-A', 'A-C', 'A-G', 'A-T', 'A-N', 'A--', ...
         'CAA', 'CAC', 'CAG', 'CAT', 'CAN', 'CA-', 'CCA', 'CCC', 'CCG', 'CCT', 'CCN', 'CC-', 'CGA', 'CGC', 'CGG', 'CGT', 'CGN', 'CG-', 'CTA', 'CTC', 'CTG', 'CTT', 'CTN', 'CT-', 'CNA', 'CNC', 'CNG', 'CNT', 'CNN', 'CN-', 'C-A', 'C-C', 'C-G', 'C-T', 'C-N', 'C--', ...
         'GAA', 'GAC', 'GAG', 'GAT', 'GAN', 'GA-', 'GCA', 'GCC', 'GCG', 'GCT', 'GCN', 'GC-', 'GGA', 'GGC', 'GGG', 'GGT', 'GGN', 'GG-', 'GTA', 'GTC', 'GTG', 'GTT', 'GTN', 'GT-', 'GNA', 'GNC', 'GNG', 'GNT', 'GNN', 'GN-', 'G-A', 'G-C', 'G-G', 'G-T', 'G-N', 'G--', ...
         'TAA', 'TAC', 'TAG', 'TAT', 'TAN', 'TA-', 'TCA', 'TCC', 'TCG', 'TCT', 'TCN', 'TC-', 'TGA', 'TGC', 'TGG', 'TGT', 'TGN', 'TG-', 'TTA', 'TTC', 'TTG', 'TTT', 'TTN', 'TT-', 'TNA', 'TNC', 'TNG', 'TNT', 'TNN', 'TN-', 'T-A', 'T-C', 'T-G', 'T-T', 'T-N', 'T--', ...
         'NAA', 'NAC', 'NAG', 'NAT', 'NAN', 'NA-', 'NCA', 'NCC', 'NCG', 'NCT', 'NCN', 'NC-', 'NGA', 'NGC', 'NGG', 'NGT', 'NGN', 'NG-', 'NTA', 'NTC', 'NTG', 'NTT', 'NTN', 'NT-', 'NNA', 'NNC', 'NNG', 'NNT', 'NNN', 'NN-', 'N-A', 'N-C', 'N-G', 'N-T', 'N-N', 'N--', ...
         '-AA', '-AC', '-AG', '-AT', '-AN', '-A-', '-CA', '-CC', '-CG', '-CT', '-CN', '-C-', '-GA', '-GC', '-GG', '-GT', '-GN', '-G-', '-TA', '-TC', '-TG', '-TT', '-TN', '-T-', '-NA', '-NC', '-NG', '-NT', '-NN', '-N-', '--A', '--C', '--G', '--T', '--N', '---'};

template_len = template_len+100;    
     
cloud1 = zeros(6,template_len);
cloud2 = zeros(36,template_len);
cloud3 = zeros(216,template_len);
 
fid = fopen(file_name);
tline= fgets(fid); 
while ischar(tline)
    
	%id = tline(1:end-1);
	pos  = sscanf(tline(2:strfind(tline,'_')-1),'%d');
    %getting seq info
    tline= fgets(fid);
    seq = tline(1:end-1);
	tline= fgets(fid);
	tline= fgets(fid);
	%quality = tline(1:end-1);
    tline = fgets(fid);

    %adding seq to the nucleotide cloud
    for i = 1:length(seq)
        temp = strcmp(nucs,seq(i));
        cloud1( temp , pos+i-1 ) = cloud1( temp , pos+i-1 ) + 1;
    end 
    %adding seq to the nucleotide pair cloud
    for i=1:(length(seq)-1)
        temp = strcmp(nucPairs,seq(i:i+1));
        cloud2( temp , pos+i-1 ) = cloud2( temp , pos+i-1 ) + 1;
    end
    %adding seq to the codon cloud
    for i=1:(length(seq)-2)
        temp = strcmp(codons,seq(i:i+2));
        cloud3( temp , pos+i-1 ) = cloud3( temp , pos+i-1 ) + 1;
    end

end

fclose(fid);

end

