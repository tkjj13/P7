clc;
close all;
clear all;

[GH GV] = gain_patch_func('patch858tot.txt',858E6,[0 pi/4])
%[GH GV] = gain_patch_func('patch24nyv2.txt',2580E6,[0 pi/4]);
%10.^(GH/10)
%10.^(GV/10)

2*pi/(3E8/858E6)
2*pi/(3E8/2580E6)


%%

sV1 = [0.01522 0.04378];
sV2 = [-0.0213 0.0436];
sV3 = [-0.02669 0.07534];

sH1 = [-0.02669 -0.00906];
sH2 = [-0.03937 0.03182];
sH3 = [-0.07578 0.02388];



n = 1000;
X = linspace(-1,1,n);
Y = X;

[x, y] = meshgrid(X,Y);

distH = sqrt((sH1(1)-x).^2+(sH1(2)-y).^2)+sqrt((sH2(1)-x).^2+...
    (sH2(2)-y).^2)+sqrt((sH3(1)-x).^2+(sH3(2)-y).^2); % horisontal pplads
distV = sqrt((sV1(1)-x).^2+(sV1(2)-y).^2)+sqrt((sV2(1)-x).^2+...
        (sV2(2)-y).^2)+sqrt((sV3(1)-x).^2+(sV3(2)-y).^2); % vertikal pplads
        
distV_vec = reshape(distV,n^2,1);
posV = find(distV_vec == min(distV_vec));

distH_vec = reshape(distH,n^2,1);
posH = find(distH_vec == min(distH_vec));


pV = X(floor(posV/n)+1)+Y(mod(posV,n))*1i;
pH = X(floor(posH/n)+1)+Y(mod(posH,n))*1i;

e0 = ((1+pV)*(1-pH))/((1+pH)*(1-pV))



%%

pplads_patch_2580_vert = -[45.207	41.497	38.173	60.751	39.732	36.045	51.803	37.628	55.86	38.293
56.899	49.82	45.098	55.157	46.772	41.415	51.41	44.79	50.539	44.448
67.641	60.209	55.091	52.979	51.059	50.511	50.905	47.452	51.021	51.222
80.334	72.054	65.765	56.403	66.967	58.658	54.883	54.538	54.848	56.864
87.852	83.383	76.121	64.371	78.76	71.413	59.642	62.597	59.145	62.113
101.535	95.524	88.822	74.713	92.376	83.113	69.565	73.48	66.442	67.517];


save pplads_meas_raw pplads_patch_858_vert pplads_patch_858_hor pplads_mono_2580_vert pplads_patch_2580_hor...
    pplads_patch_2580_vert pplads_mono_858_vert pplads_mono_858_hor pplads_mono_2580_vert pplads_mono_2580_hor pplads_demo_vert pplads_demo_hor;




