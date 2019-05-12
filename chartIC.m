function [ mT2,sd,Ttr,Itr] = chartIC( T2tr,varargin)
% build in control chart

[del,ifIncludeNAN] = process_options(varargin,'del',1,'ifIncludeNAN',1);

T2trd= T2tr(:,(del+1):end);

if ifIncludeNAN
    mT2 = nanmean(T2trd,2);
    sd = nanstd(T2trd,0,2);
else
    mT2 = mean(T2trd,2);
    sd = std(T2trd,0,2);
end

T2c = bsxfun(@minus,T2trd,mT2);
T2tn = bsxfun(@rdivide,T2c,sd);
[Ttr,Itr]=(max(T2tn,[],1));


end

