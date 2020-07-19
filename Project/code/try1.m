bps=4;
M=2^bps;
nBits=1e3*bps;
ebno=10;


symMap       = [11,10,14,15,9,8,12,13,1,0,4,5,3,2,6,7];
bitTable     = de2bi(symMap,bps,'left-msb');

sym          = qammod(symMap(1:M),M,symMap,'UnitAveragePower',false,'PlotConstellation',true);
sym          = reshape(sym,[M,1]);
averagePower = sum(abs(sym).^2)/16;


Rayleigh = comm.RayleighChannel('PathGainsOutputPort',true,'RandomStream','Global stream','NormalizePathGains',true);
awgnChan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)','SNR',10,'SignalPower',averagePower);
sphDec   = comm.SphereDecoder('Constellation',sym,'BitTable',bitTable,'DecisionType','Hard');
berRate  = comm.ErrorRate;


data            = randi([0 1],nBits,1);
modData         = qammod(data,M,symMap,'InputType','bit','UnitAveragePower',false,'PlotConstellation',true);
[y,pathGains]   = step(Rayleigh,modData);
rxsig           = step(awgnChan,y);
decodedData     = step(sphDec,rxsig,pathGains); 
dataOut         = double(decodedData(:));

errorStats = step(berRate,data,dataOut);
errorStats(1:2)