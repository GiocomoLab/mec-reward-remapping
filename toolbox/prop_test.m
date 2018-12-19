function [z,p] = prop_test(n1,n2,total1,total2)
 
p2 = n1/total1;
p1 = n2/total2;
p = (n1 + n2)/(total1+total2);
z = abs((p1-p2)/sqrt(p*(1-p)*((1/total1)+(1/total2))));
p = 2*(1-normcdf(z,0,1));
 
return