function rxsig = awgn(y1,snr,averagePower)
    %variance of complex noise
    N0 = averagePower/snr;
    
    %std. deviation of real and i component
    sigma = sqrt(N0/2);
    
    %real + imaginary part
    n1 = normrnd(0,sigma,[length(y1),1]);
    n2 = normrnd(0,sigma,[length(y1),1]);
    
    %complex noise
    n = n1 + 1j*n2;
    
    %received signal at rx
    rxsig = y1 + n;

end
