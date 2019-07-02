function A=cgl(N,c1,c3,M)
% A = cgl(N,omega,tol)
% Pseudo-spectral solution of complex Ginsburg-Landau equation
% N = number of grid points in both dimensions
% c1,c3 = equation coefficients
% M = number of time steps

L = 128*pi;
T = 10000;
dt = T/M;

% initialize A with random data
A = 3*rand(N,N)-1.5+i*(3*rand(N,N)-1.5);

figure(1)
subplot(2,2,1)
contourf(abs(A));
title(sprintf('time = %f',0));
axis equal
subplot(2,2,2)
contourf(atan2(imag(A),real(A)),[-pi:pi/20:pi]);
axis equal
subplot(2,2,3)
surf(abs(A));
subplot(2,2,4)
surf(atan2(imag(A),real(A)));
axis([0,N,0,N,-pi,pi]);
rng=[];
ptime = tic;

for n=1:T/dt
    % compute RHS
    A1 = A+dt/4*(A+(1+i*c1)*del2A(A)*(2*pi/L)^2-(1-i*c3)*abs(A).^2.*A);
    A2 = A+dt/3*(A1+(1+i*c1)*del2A(A1)*(2*pi/L)^2-(1-i*c3)*abs(A1).^2.*A1);
    A1 = A+dt/2*(A2+(1+i*c1)*del2A(A2)*(2*pi/L)^2-(1-i*c3)*abs(A2).^2.*A2);
    A = A+dt*(A1+(1+i*c1)*del2A(A1)*(2*pi/L)^2-(1-i*c3)*abs(A1).^2.*A1);
    
    if toc(ptime) >= 1
        ptime = tic;
        subplot(2,2,1)
%        contourf(abs(A),[0:0.02:1]);
        contourf(abs(A));
        title(sprintf('time = %f',n*dt));
        axis equal
        subplot(2,2,2)
        contourf(atan2(imag(A),real(A)),[-pi:pi/20:pi]);
        axis equal
        subplot(2,2,3)
        surf(abs(A));
        subplot(2,2,4)
        surf(atan2(imag(A),real(A)));
        axis([0,N,0,N,-pi,pi]);
        drawnow
    end
end
end

function d2A = del2A(A)
N = size(A,1)/2;
d2A = fftshift(fft2(A));
d2Ax = -(ones(2*N,1)*[-N:N-1].^2).*d2A;
d2Ay = -d2A.*(([-N:N-1]'*ones(1,2*N)).^2);
d2A = ifft2(fftshift(d2Ax+d2Ay));
end