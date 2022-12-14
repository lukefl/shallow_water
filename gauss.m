function y=gauss(mean,x, sig)
n=size(x,2);
y=zeros(1,n);
y=(1/sqrt(2*pi*sig)).*exp(-0.5.*((x-mean)/sig).^2);