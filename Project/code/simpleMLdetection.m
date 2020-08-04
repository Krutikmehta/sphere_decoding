function decodedData = simpleMLdetection(rxsig,H,M,symMap)

%constellation
bitTable     = de2bi(symMap,4,'left-msb');

%symbolset modulation and reshaping
Re = -(sqrt(M)-1):2:(sqrt(M)-1);
Im = 1j*Re;
%signal space
sym = repmat(Re,4,1);
sym = sym(:) + repmat(Im,1,4)';

len = length(rxsig);
diff = zeros([M,M]);
decodedData = zeros(len,1);
for k =1:len/2
    
    for i =1:M
        for j = 1:M
            diff(i,j) = vecnorm(rxsig(2*k-1:2*k)-H*[sym(i);sym(j)]);
        end
    end
    minMatrix = min(diff(:));
    [row,col] = find(diff==minMatrix);
    decodedData(2*k-1:2*k) = [sym(row);sym(col)];
end

end