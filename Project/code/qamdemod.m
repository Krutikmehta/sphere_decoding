function demodData = qamdemod(data,M,symMap)

Re = -(sqrt(M)-1):2:(sqrt(M)-1);
Im = 1j*Re;

%signal space
sigspace = repmat(Re,4,1);
sigspace = sigspace(:) + repmat(Im,1,4)';

symbs = zeros(size(data,1),1);

for i = 1:size(data,1)
    symbs(i) = symMap(double(find(sigspace==data(i))));
    
end

demodData = de2bi(symbs,sqrt(M));
demodData = flip(demodData,2);

end