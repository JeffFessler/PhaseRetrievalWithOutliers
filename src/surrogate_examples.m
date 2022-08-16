% surrogate_examples.m -- This code generates the figures comparing the
% 1-norm data fit term against its convex majorizer, shown in Fig. 1 in the
% paper.
%
% This code is subject to copyright and the license set forth in
% LICENSE.TXT. If you did not receive a copy of LICENSE.TXT with this
% software, or have other questions about the code, please contact Daniel
% Weller (University of Virginia) at d.s.weller@ieee.org.

tmax = 2;
y = 1;
q = 2;
s1 = -0.5;
s2 = 1.25;

hfun = @(t) abs(y-abs(t).^q);
phifun1 = @(t,s) max(abs(t).^q-y,y+(q-1).*abs(s).^q-q.*abs(s).^(q-1).*real(t.*exp(complex(0,-angle(s)))));
phifun = @(t,s) phifun1(t,min(abs(s),nthroot(y,q)).*exp(complex(0,angle(s))));

ts = linspace(-tmax,tmax,101).';
hs = hfun(ts);
phis1 = phifun(ts,s1);
phis2 = phifun(ts,s2);
phis2b = phifun1(ts,s2);

figure;
subplot(1,2,1);
plot(ts,[hs,phis1],'-',s1,hfun(s1),'o');

subplot(1,2,2);
plot(ts,[hs,phis2,phis2b],'-',s2,hfun(s2),'o');
