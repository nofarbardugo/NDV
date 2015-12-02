function [tcon, full_tcon] = tcloud2con(tcloud)
% gets consensus from cloud of all contings.
% tcon = resulting consensus without gaps
% full_tcon = resulting consensus with gaps

tcon='-';
for i=1:size(tcloud,2)
    tcon(i)='-';
    for j=1:size(tcloud,1)
        if tcloud(j,i)~='-'
            tcon(i) = tcloud(j,i);
        end
    end
end

full_tcon = tcon;
tcon(tcon=='-')=[];

end