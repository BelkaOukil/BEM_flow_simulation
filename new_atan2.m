function res = new_atan2(y,x)
res = atan2(y,x);
if x<0 && y<0
    res = res + 2*pi;
end