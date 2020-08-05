function decodedData = untitled(rxsig,H,M,symMap)


%symbolset modulation and reshaping
Re = -(sqrt(M)-1):2:(sqrt(M)-1);
Im = 1j*Re;
%signal space
sym = repmat(Re,sqrt(M),1);
sym = sym(:) + repmat(Im,1,sqrt(M))';
%sym          = reshape(sym,[M,1]);
averagePower = sum(abs(sym).^2)/16;
len= length(rxsig);


rxsig  = reshape(rxsig,[2,len/2]);

a = repmat(sym,[M,1]);
b = repmat(sym,[1,M])';
signalspace = [b(:),a]';
recievespace = H*signalspace;
decodedData = zeros(2,len/2);

for i =1:len/2
    diff = vecnorm(abs(recievespace - rxsig(:,i)));
    [M1,I] = min(diff,[],2);
    %decodedData(1,i) = signalspace(1,I);
    %decodedData(2,i) = signalspace(2,I);
    decodedData(:,i) = signalspace(:,I);
end

decodedData = decodedData(:);

end
