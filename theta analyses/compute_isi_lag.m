function [Nm,m] = compute_isi_lag(cellTS)

isi_lag = repmat(cellTS,2001,1)'-lagmatrix(cellTS,(-1000:1000));
isi_lag = isi_lag/1e3;
maxLag = 500;
isi_lag(abs(isi_lag) > maxLag) = [];
isi_lag(isi_lag == 0) = [];
hist_bin_size = 10; % 3 milliseconds
[Nm,m] = hist(isi_lag,maxLag*2/hist_bin_size);
Nm = (Nm - min(Nm))./range(Nm);

return