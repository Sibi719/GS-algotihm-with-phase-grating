clc;
clear all;
M=512;
N=512;
z=0.02;
Tp=(635*10^-9);
d1=(10*10^-3)/M;
it=40;
fm=1/Tp;
lambda=635*10^-9;
noise=randn(N).*sqrt(10^-4);
wn=1;
k=(2*pi)/lambda;
xm= (-M/2:M/2-1)*d1;
ym= (-M/2:M/2-1)*d1;
[Xm,Ym]=meshgrid(xm,ym);
xn= (-N/2:N/2-1)*d1;
yn= (-N/2:N/2-1)*d1;
[Xn,Yn]=meshgrid(xn,yn);

S=zeros(N);
R=(M/2)*d1;
A=abs(Xn)<=R&abs(Yn)<=R;
S(A)=1;

figure
imagesc(xn*10^3,yn*10^3,S);
colormap(gray);
xlabel('x(mm)');
ylabel('y(mm)');
colormap(gray);
colorbar
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"1.png");

T=imread("Lena.png");
T=rgb2gray(T);
T=double(T);
T=rescale(T,0.3,0.92);
Nr=size(T,1);
Nc=size(T,2);
Dr=(N-Nr)/2;
Dc=(N-Nr)/2;
T = padarray(T,[Dr,Dc],0,'both');
figure
imagesc(xn*10^3,yn*10^3,T);
xlabel('x(mm)');
ylabel('y(mm)');
colormap(gray);
colorbar
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"2.png");

P=imread("dog.jpg");
P=rgb2gray(P);
P=double(P);
P=rescale(P,-0.4*pi,0.3*pi);
Nr=size(P,1);
Nc=size(P,2);
Dr=(N-Nr)/2;
Dc=(N-Nr)/2;
P = padarray(P,[Dr,Dc],0,'both');
figure
imagesc(xn*10^3,yn*10^3,P);
xlabel('x(mm)');
ylabel('y(mm)');
colormap(gray);
colorbar
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"3.png");



Uo=ones(M);
a=0.65;
n=0.5;
theta1=0;
theta2=90;


A= pi*0.5*(square(2.*pi.*(fm*cosd(theta1).*Xm+fm*sind(theta1).*Ym))+1);
A = padarray(A,[Dr,Dc],0,'both');

B= pi*0.5*(square(2.*pi.*(fm*cosd(theta2).*Xm+fm*sind(theta2).*Ym))+1);
B = padarray(B,[Dr,Dc],0,'both');

MA=exp(1j.*A);
MB=exp(1j.*B);

figure
imagesc(xn*10^3,yn*10^3,A);
xlabel('x(mm)');
ylabel('y(mm)');
colorbar
colormap(gray)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"4.png");

figure
imagesc(xn*10^3,yn*10^3,B);
colorbar
colormap(gray)
xlabel('x(mm)');
ylabel('y(mm)');
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"5.png");



U1A=sqrt(T).*exp(1j.*P).*MA;
U1B=sqrt(T).*exp(1j.*P).*MB;

[Uo,Io]=FresProp(d1,0,lambda,N,sqrt(T).*exp(1j.*P));
[U1,I1]=FresProp(d1,z,lambda,N,sqrt(T).*exp(1j.*P));
[U2A,I2A]=FresProp(d1,z,lambda,N,U1A);%forward_Fresnel(k,z,lambda,Xn,Yn,U1A);
[U2B,I2B]=FresProp(d1,z,lambda,N,U1B);%forward_Fresnel(k,z,lambda,Xn,Yn,U1B);

Io=Io+ (wn.*noise);
I1= I1+ (wn.*noise);
I2A = I2A + (wn.*noise);
I2B = I2B + (wn.*noise);


figure
imagesc(xn*10^3,yn*10^3,Io);
colorbar
colormap(gray)
xlabel('x(mm)');
ylabel('y(mm)');
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];


figure
imagesc(xn*10^3,yn*10^3,I1);
colorbar
colormap(gray)
xlabel('x(mm)');
ylabel('y(mm)');
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];


figure
imagesc(xn*10^3,yn*10^3,I2A);
colorbar
colormap(gray)
xlabel('x(mm)');
ylabel('y(mm)');
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

figure
imagesc(xn*10^3,yn*10^3,I2B);
colorbar
colormap(gray)
xlabel('x(mm)');
ylabel('y(mm)');
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

AA=sqrt(Io);
BB=sqrt(I1);
CC=sqrt(I2A);
DD=sqrt(I2B);


%%GS algothm


U1_est_abs = AA.*S;
U1_est_phase= rand(N).*S;

U1_est =U1_est_abs.*exp(1j.*U1_est_phase);

for i=1:it

 [U2_est,~]=FresProp(d1,z,lambda,N,U1_est);
 U2_est= BB.*exp(1j.*angle(U2_est));
 [U1_est,~] = FresProp(d1,-z,lambda,N,U2_est);
 U1_est= AA.*exp(1j.*angle(U1_est));
 
  U1_est_A=U1_est.*MA;
 [U2_est_A,~]=FresProp(d1,z,lambda,N,U1_est_A);
 U2_est_A= CC.*exp(1j.*angle(U2_est_A));
 [U1_est_A,~] = FresProp(d1,-z,lambda,N,U2_est_A);
 U1_est=U1_est_A./MA;
 U1_est = AA.*exp(1j.*(angle(U1_est)));
 
 U1_est_B=U1_est.*MB;
 [U2_est_B,~]=FresProp(d1,z,lambda,N,U1_est_B);
 U2_est_B= DD.*exp(1j.*angle(U2_est_B));
 [U1_est_B,~] = FresProp(d1,-z,lambda,N,U2_est_B);
 U1_est=U1_est_B./MB;
 U1_est = AA.*exp(1j.*(angle(U1_est)));
 
  [Uee,Iee]=FresProp(d1,z,lambda,N,U1_est);
 
  err=((Iee - I1)./I1).^2;
  error(i)= (sum(sum(err)*d1*d1));
end

figure
imagesc(xn*10^3,yn*10^3,(abs(U1_est)).^2);
xlabel("x(mm)");
ylabel("y(mm)");
colormap(gray)
colorbar
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"12.png");


figure
imagesc(xn*10^3,yn*10^3,(angle(U1_est)));
xlabel("x(mm)");
ylabel("y(mm)");
colormap(gray)
colorbar
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"13.png");

i=1:it;
% figure
% plot(i,error);
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"14.png");

figure
plot(i,log(error));
xlabel("Iterations");
ylabel("Error");
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf,"15.png");
% 
% figure
% imagesc(xn*10^3,yn*10^3,abs(U2A));
% colormap(gray);
% xlabel('x(mm)');
% ylabel('y(mm)');
% colormap(gray);
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"16.png");
% figure
% imagesc(xn*10^3,yn*10^3,abs(U2B));
% colormap(gray);
% xlabel('x(mm)');
% ylabel('y(mm)');
% colormap(gray);
% colorbar
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% saveas(gcf,"17.png");

function [U,I]=forward_Fresnel(k,z,lambda,Xn,Yn,Uin)
 U=((exp(1j.*k.*z)/(1j*lambda*z)).*exp(((1j*k)/(2*z)).* (Xn.^2+Yn.^2)).*fftshift(fft2(Uin.*exp(((1j*k)/(2*z)).*(Xn.^2+Yn.^2)))));
 I=abs(U).^2;
end

function [U,I]=back_Fresnel(k,z,lambda,Xn,Yn,Uin)
 U=((1j*lambda*z)./exp(1j.*k.*z)).*exp(((-1j*k)/(2*z)).*(Xn.^2+Yn.^2)).*ifft2(ifftshift(Uin.*exp(((-1j*k)/(2*z)).*(Xn.^2+Yn.^2))));
 I=abs(U).^2;
end

function [Ri,R]=FresProp(dpix,d,lambda,Hsize,Hcrop)
 z=d;

%Spatial frequencies
Xsize = Hsize*dpix; %Hsize is the number of pixel, dpix is the length of one pixel, Xsize is the total lenght of the image. 
du = 1/(Xsize);% the resolution of fourier frequency coordinates
%Nyquist cut-off for Sampling Hologram
umax = 1/(2*dpix); %define the k space 
u = -umax:du:umax-du;
[U,V]=meshgrid(u,u);
clear  u V  du;
 
%Evanescent cut-off 
uev = 1/lambda; %???????
 
%Nyquist cut-off for Fresnel Propagation Kernel
unp = uev*(Xsize/(2*abs(z)));
clear Xsize;
 
%Circular window
A = U.^2+(U').^2;
clear U;
if uev>=unp
    ucut = unp;
end
if unp>uev
    ucut = uev;
end
W= sqrt(A);
W = (W<=ucut); 
% disp(['Cutoff =',num2str(ucut),' Evansecent Cutoff =',num2str(uev),...
%' Nyquist Cutoff =', num2str(unp),'u max =',num2str(umax)])
clear ucut uev unp
 
%Fresnel kernel: paraxial approximation
H = exp((-1i*pi*lambda* z).*(A));
clear A;
 
%Truncate kernel
H = W.*H;
clear W;
 
%Hologram Spectrum
Htemp = fft2(Hcrop);
HH = fftshift(Htemp);
clear Htemp;
 
%Propagate field
RR = HH.*H;
clear H HH;
RR =ifftshift(RR);
Ri = ifft2(RR); 
R=abs(Ri).^2;

end
