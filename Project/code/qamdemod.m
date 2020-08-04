function demodData = qamdemod_idp(data,M,symMap)
%M = 16;
%symMap       = [11 10 14 15 9 8 12 13 1 0 4 5 3 2 6 7];
%data            = [1.000000000000000 - 3.000000000000000i;3.000000000000000 + 3.000000000000000i;3.000000000000000 - 1.000000000000000i;1.000000000000000 + 3.000000000000000i];

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