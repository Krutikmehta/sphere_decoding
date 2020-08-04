bps=4; %bits per signal
M=2^bps;
nBits=1e3*bps;
ebno=10; %SNR

%constellation
symMap       = [11 10 14 15 9 8 12 13 1 0 4 5 3 2 6 7];
bitTable     = de2bi(symMap,bps,'left-msb');

%symbolset modulation and reshaping
Re = -(sqrt(M)-1):2:(sqrt(M)-1);
Im = 1j*Re;
%signal space
sym = repmat(Re,4,1);
sym = sym(:) + repmat(Im,1,4)';
%sym          = reshape(sym,[M,1]);
averagePower = sum(abs(sym).^2)/16;


%random input datastream
data            = randi([0 1],nBits,1);
%modulation
modData = qammod(data,M,symMap);


h11_r = normrnd(0,(1/sqrt(2)));
h11_i = normrnd(0,(1/sqrt(2)));
h12_r = normrnd(0,(1/sqrt(2)));
h12_i = normrnd(0,(1/sqrt(2)));
h21_r = normrnd(0,(1/sqrt(2)));
h21_i = normrnd(0,(1/sqrt(2)));
h22_r = normrnd(0,(1/sqrt(2)));
h22_i = normrnd(0,(1/sqrt(2)));
h11 = h11_r + 1j*h11_i;
h12 = h12_r + 1j*h12_i;
h21 = h21_r + 1j*h21_i;
h22 = h22_r + 1j*h22_i;

H = [[h11,h12];[h21,h22]]; %matrix H
len=nBits/bps; %Number of signals
x = reshape(modData,[2,len/2]);
y1 = H*x;
y1 = y1(:);

rxsig           = awgn(y1,ebno,averagePower);
decodedData = zeros(len,1);

for i =1:len
    diff = abs(sym - rxsig(i))
    [M1,I] = min(diff);
    decodedData(i) = sym(I);
    
end


%demodulation
rxData = qamdemod(decodedData,M,symMap)';
rxData = rxData(:);

errorStats = ber(data,rxData);