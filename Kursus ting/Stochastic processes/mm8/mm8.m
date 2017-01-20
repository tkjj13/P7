clc
clear all
close all


randn(1);
h = 500;
U = randn(1);
Uold = randn(1);

for k = 1:10
  for n = 1:h 
    x(k,n) = U + Uold;
    Uold = U;
    U = randn(1);
  end


for n = 1:h 
    y(k,n) = U - Uold;
    Uold = U;
    U = randn(1);
end
end

for s = 1:k
%%biased
for n = 1:h-1
    g = 0;
    for l = 1:h-n
        g = x(s,l)*x(s,l+n-1)+g;
    end    
acfx(s,n) = 1/h*g;
end


%%unbiased
for n = 1:h-1
    g = 0;
    for l = 1:h-n
        g = x(s,l)*x(s,l+n-1)+g;
    end    
acfx1(s,n) = (1/(h-n))*g;
end
end

acfxres = sum((1/k)*(acfx));
acfx1res = sum((1/k)*(acfx1));





stem(acfx1(1,:));
hold on
stem(acfx(1,:))
figure;
return
autocorr(x);
autocorr(y);
return
plot(abs(fft(x)));
plot(abs(fft(y)));



