bps=4;
M=2^bps;
nBits=1e3*bps;
ebno=10;


symMap       = [11,10,14,15,9,8,12,13,1,0,4,5,3,2,6,7];
bitTable     = de2bi(symMap,bps,'left-msb');

sym          = qammod(symMap(1:M),M,symMap,'UnitAveragePower',false,'PlotConstellation',false);
sym          = reshape(sym,[M,1]);
averagePower = sum(abs(sym).^2)/16;


Rayleigh = comm.RayleighChannel('PathGainsOutputPort',true,'RandomStream','Global stream','NormalizePathGains',true);
awgnChan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)','SNR',10,'SignalPower',averagePower);
sphDec   = comm.SphereDecoder('Constellation',sym,'BitTable',bitTable,'DecisionType','Hard');
berRate  = comm.ErrorRate;


data            = randi([0 1],nBits,1);
modData         = qammod(data,M,symMap,'InputType','bit','UnitAveragePower',false,'PlotConstellation',false);
%[y,pathGains]   = step(Rayleigh,modData);
for 
rxsig           = step(awgnChan,y);

rxsig_real      = real(rxsig);
rxsig_img       = imag(rxsig);

len=nBits/bps;
for j=1:len
    rxsig_2(2*j-1)=rxsig_real(j);
    rxsig_2(2*j)=rxsig_img(j);
end

rxsig_2=reshape(rxsig_2,[2*len,1]);

h11_r = normrnd(0,(1/sqrt(2)));
h11_i = normrnd(0,(1/sqrt(2)));
h12_r = normrnd(0,(1/sqrt(2)));
h12_i = normrnd(0,(1/sqrt(2)));
h21_r = normrnd(0,(1/sqrt(2)));
h21_i = normrnd(0,(1/sqrt(2)));
h22_r = normrnd(0,(1/sqrt(2)));
h22_i = normrnd(0,(1/sqrt(2)));

B = [[h11_r,-h11_i,h12_r,-h12_i];[h11_i,h11_r,h12_i,h12_r];[h21_r,-h21_i,h22_r,-h22_i];[h21_i,h21_r,h22_i,h22_r]];
A = zeros(4,1)+3;

y_c = rxsig_2 + B*A;
[Q,R] = qr(2*B);


decodedData     = step(sphDec,rxsig,pathGains); 
dataOut         = double(decodedData(:));

errorStats = step(berRate,data,dataOut);
errorStats(1:2)