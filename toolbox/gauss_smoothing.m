function x=gauss_smoothing(x,sigma)
if sigma>0
    % define gaussian filter
    lw=ceil(3*sigma);
    wx=-lw:lw;
    gw=exp(-wx.^2/(2*sigma^2)); gw=gw/sum(gw);

    x=conv(x,gw,'same');
end

