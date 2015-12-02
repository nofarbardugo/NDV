
function [tcloud, tcon] = deNovo_check(long_seq, seq_ln, seq_fit_per,newcastle_con)
% combines deNovo built sequenes to one contig.
% long_seq = deNovo resulting sequences
% seq_len = minimal sequence length to use
% seq_fit_per = minimal % of fit to known sequnce to use
% newcastle_con = known viral sequence for compare
% tcloud = a cloud made of all related sequences combined
% tcon = consensus sequence of sequnces combined

clear a b c tcloud pos

% remove short seuqnces
for i=length(long_seq):-1:1
    if length(long_seq{i}) < seq_ln
        long_seq(i) = [];
    end
end

%making basic cloud
tcloud ='-';
for i=1:length(long_seq)
    for j=1:length(newcastle_con)
        tcloud(i,j)='-';
    end
end

pos=[];
bad=[];
for i=1:length(long_seq)
    
    % finding right alignment (forward or reverse)
    [a1 b1 c1] = swalign(seqrcomplement(long_seq{i}),newcastle_con);
    [a2 b2 c2] = swalign(long_seq{i},newcastle_con);
    
    if a1>a2
        tc=c1(2);
        tb=b1;
    else
        tc=c2(2);
        tb=b2;
    end
    
    %checking fit to known sequence
    if sum(tb(2,:)~='|')/length(tb) > seq_fit_per
        bad(end+1)=i;
    end
    
    %adding to cloud
    pos(1,i) = tc;
    pos(2,i) = tc + length(tb)-1;
    
    tcloud(i,pos(1,i):pos(2,i)) = tb(1,:);
    
    le(i)=length(tb);
    
end
% removing bad sequences
tcloud(bad,:)=[];
long_seq(bad)=[];
pos(:,bad)=[];
%%
% finding consensus
tcon='-';
for i=1:size(tcloud,2)
    tcon(i)='-';
    for j=1:size(tcloud,1)
        if tcloud(j,i)~='-'
            tcon(i) = tcloud(j,i);
        end
    end
end
tcon(tcon=='-')=[];

end
