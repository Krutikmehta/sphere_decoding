function modData = qammod(data,M,symMap)

data = reshape(data,[sqrt(M), size(data,1)/sqrt(M)]);
data = data';
symbs = bi2de(data,'left-msb');

Re = -(sqrt(M)-1):2:(sqrt(M)-1);
Im = 1j*Re;

%signal space
sigspace = repmat(Re,4,1);
sigspace = sigspace(:) + repmat(Im,1,4)';

modData = zeros(size(symbs,1),1);
%symMap = symMap + 1;

for i = 1:size(symbs,1)
    modData(i) = double(symMap==symbs(i))*sigspace;
    
end

end