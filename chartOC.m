function [ Ttr,Itr] = chartOC( T2tr,mT2,sd)
% build in control chart
%[del] = process_options(varargin,'del',1);


T2c = bsxfun(@minus,T2tr,mT2);
T2tn = bsxfun(@rdivide,T2c,sd);
[Ttr,Itr]=(max(T2tn,[],1));


end

