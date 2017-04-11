function flag = isFalsePositive(muRxyz)

flag = false;
d = muRxyz - [0 0 -40]';  % to remove temporal pole
r1 = sqrt(sum(d.*d));
if r1<30, flag=true; end;

d = muRxyz - [0 -74 70]';  % to remove outside brain (above SMA)
r2 = sqrt(sum(d.*d));
if r2<10, flag=true; end;

d = muRxyz - [-10 -105 10]';  % to remove outside brain (posterior visual)
r3 = sqrt(sum(d.*d));
if r3<10, flag=true; end;

d = muRxyz - [10 -105 10]';  % to remove outside brain (posterior visual)
r4 = sqrt(sum(d.*d));
if r4<10, flag=true; end;

% d = muRxyz - [-22 -42 12]';  % to remove left ventricle
% r2 = sqrt(sum(d.*d));
% if r2<8, flag=true; end;
% 
% d = muRxyz - [22 -42 12]';  % to remove right ventricle
% r3 = sqrt(sum(d.*d));
% if r3<8, flag=true; end;
