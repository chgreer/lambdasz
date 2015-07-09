function w = tophat_window(x)
%
%function w = tophat_window(x)
%
%Represents the fourier-transform of a real-space 3D tophat.
%
%        w = ( 3.0.*(sin(x)-x.*cos(x))./(x.*x.*x) );
%

w = ( 3.0.*(sin(x)-x.*cos(x))./(x.*x.*x) );

return

