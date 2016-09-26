function hpx = highPass(x)
% hps = highPass(x)
%
% Simple high-pass filtering routine. 
    org_size = size(x);
    qmf = MakeONFilter('Daubechies',8);
    x = padarray(x,[8, 8],'symmetric');
    spectrum = mdwt(x,qmf,2);
    region = ceil(size(x)/4);
    spectrum(1:region(1), 1:region(2)) = 0;
    hpx = midwt(spectrum, qmf, 2);
    hpx = imcropmiddle(hpx, org_size);
end
