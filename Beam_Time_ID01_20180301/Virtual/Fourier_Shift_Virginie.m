clear;
close all;
N1 = 400;
N2 = 400;
phi = 2*pi;%2*pi/4;
[r2,r1] = meshgrid((1:1:N1)-N1/2, (1:1:N2)-N2/2);
 
object = (abs(r2 -15)<=20).*(abs(r1)<=10);
figure(1);
imagesc(object); title('object');
Fourier = fftshift(fftn(fftshift(object)));
Fourier_shift = fftshift(fftn(fftshift(object.*exp(1i*phi*r2))));
 
figure(2);
subplot(2,2,1);
imagesc(log10(abs(Fourier))); title('Fourier spectrum');
subplot(2,2,2);
imagesc(angle(Fourier));
subplot(2,2,3);
imagesc(log10(abs(Fourier_shift))); title('Shifted Fourier spectrum ');
subplot(2,2,4);
imagesc(angle(Fourier_shift));
 
 
back_Fourier = fftshift(ifftn(fftshift(Fourier)));
back_Fourier_shift = fftshift(ifftn(fftshift(Fourier_shift)));
 
figure(3);
subplot(2,2,1);
imagesc(log10(abs(back_Fourier))); caxis([-8 0]); colorbar; title('Back Fourier spectrum');
subplot(2,2,2);
imagesc(angle(back_Fourier));
subplot(2,2,3);
imagesc(log10(abs(back_Fourier_shift))); caxis([-8 0]); colorbar; title('Back Shifted Fourier spectrum ');
subplot(2,2,4);
imagesc(angle(back_Fourier_shift));
