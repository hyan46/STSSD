function [T2,Snow,thetanow,t,Itr,defect,Tte,gammaauto] = ewmamonit( Y,B,Bs,lambda,allgamma,varargin)
% parameters
lambdat = lambda(end);
lambdaxy = lambda(1:end-1);
% lambdat1 = 1/(1+lambdat);
lambdat1 = lambdat;
% process additional parameters
[isewma,maxIter,type,L,mT2,sd,issave,initial] = process_options(varargin,'isewma',0,'maxIter',1,'type','s','L',0,'mT2',0,'sd',1,'issave',0,'initial',0);

sizey = size(Y);
ndim = length(sizey);


if length(lambdaxy) == 1
    lambdaxy = lambdaxy*ones(ndim-1,-3);
end


% define additional functions
vec = @(x) x(:);
gammalength = length(allgamma);
constructD =@(n) diff(eye(n),1);

if iscell(Y)
    nT = length(Y);
    chooseI = @(Y,i) Y{i};
else
    nDim = length(size(Y));
    nT = size(Y,nDim);
end

if issave
    
    thetaall = zeros(size(Y,1),size(Y,2),nT);
    Sall = zeros(size(Y,1),size(Y,2),nT);
end

if isempty(B)
    for idim = 1:(ndim-1)
        B{idim} = eye(size(Y,idim));
    end
end


if isempty(Bs)
    LL = 2;
    isOrth = 1;
    BetaS = Y*0;
elseif length(Bs) == 2
    LL = 2*norm(Bs{1})^2*norm(Bs{2})^2;
    isOrth = 2;
    X = zeros(size(Bs{1},2),size(Bs{2},2));
    BetaS = zeros(size(Bs{1},2),size(Bs{2},2));

end

% preprocessing
D = cell(ndim-1,1);
K = cell(ndim-1,1);
H = cell(ndim-1,1);

for idim = 1:(ndim-1)
    if isempty(B{idim})
        B{idim} = eye(size(Y,idim));
    end
    D{idim} = constructD(size(B{idim},2));
    H{idim} = B{idim}/(B{idim}'*B{idim} + lambdaxy(idim) * (D{idim}'*D{idim}))*B{idim}';
end

if iscell(Y)
    nT = length(Y);
    chooseI = @(Y,i) Y{i};
else
    nDim = length(size(Y));
    nT = size(Y,nDim);
    if nDim == 2
        chooseI = @(Y,i) Y(:,i);
    elseif nDim == 3
        chooseI = @(Y,i) Y(:,:,i);
    elseif nDim == 4
        chooseI = @(Y,i) Y(:,:,:,i);
    else
        error('nDim not surpported,please added in here')
    end
end



T2 = zeros(gammalength,nT);
Tte = zeros(1,nT);
Itr = zeros(1,nT);
defect=0;
tnew = 1;
                   

for t = 1:nT
    tic
    dall = cell(gammalength,1);
    thetai = cell(gammalength,1);
    e = zeros(gammalength,1);
    for i = 1:gammalength
        %
        y = chooseI(Y,t);
        if t==1
            if ndim == 2
                thetanow = H{1}*(chooseI(Y,1));
            elseif ndim == 3
                thetanow = H{1}*(chooseI(Y,1))*H{2};
            elseif ndim >= 4
                thetanow = double(ttm(tensor(chooseI(Y,1)),H));
            end
            %if issave
            %    thetaall{i,1} = thetanow;
            %end
        else
            Snow = 0;
            Snowold = 1;
            iiter = 0;
            thetaold = thetanow;
            while iiter <maxIter
                iiter = iiter +1;
                Snowold = Snow;
                BetaSold = BetaS;
                told = tnew;
                if ~isewma
                    if ndim == 2
                        yhat = H{1}*(y-Snow);
                    elseif ndim == 3
                        yhat = H{1}*(y-Snow)*H{2};
                    elseif ndim >= 4
                        yhat = double(ttm(tensor(y-Snow),H));
                    end
                    thetanow = lambdat1*yhat+(1-lambdat1)*thetaold;
                else
                    thetanow = lambdat1*(y-Snow)+(1-lambdat1)*thetaold;
                end
                if isempty(Bs)
                    Snow = wthresh(y - thetanow,type,allgamma(i));                
                else
                    BetaSe = X + 2/LL*Bs{1}'*(y -Bs{1}*X*Bs{2}' - thetanow)*Bs{2};
                    BetaS = wthresh(BetaSe,type,allgamma(i)/LL); % sign(BetaSe) .* plus0(abs(BetaSe)- gamma/L);
                    Snow = Bs{1} *BetaS* Bs{2}';
                    tnew = (1+sqrt(1+4*told^2))/2;
                    if iiter==1
                        X = BetaS;
                    else
                        X = BetaS+(told-1)/tnew*(BetaS-BetaSold);
                    end
                end
            end
%            if issave
%                thetaall{i,t} = thetanow;
%           end
        end
        if isempty(Bs)
            Ye = (y - thetanow);
            maxYe = max(abs(Ye(:)));
            gammaauto = 2*graythresh(abs(Ye)/maxYe)*maxYe;
            d = wthresh(Ye,type,allgamma(i));
            
        else
            BetaSe = BetaS + 2/LL*Bs{1}'*(y -Bs{1}*BetaS*Bs{2}' - thetanow)*Bs{2};
            BetaS = wthresh(BetaSe,type,allgamma(i)/LL); % sign(BetaSe) .* plus0(abs(BetaSe)- gamma/L);
            Snow = Bs{1} *BetaS* Bs{2}';
            d = Snow; 
        end
        
        T2(i,t) = (sum(vec(d.*(y - thetanow))))^2/(sum(vec(d.*d)));
        %if issave
        %    Sall{i,t} = d;
        %end
        dall{i} = d;
        thetai{i} = thetanow;
        e(i) = sum(sum((y-d-thetanow).^2)) + 2*sum(sum(d~=0));
    end

    if L % Check for Control Limit
        [L1,Itr(t)] = chartOC( T2(:,t),mT2,sd);
        Tte(t) = L1;
        if L1 > L &&t>initial
            defect = dall{Itr(t)};
            break;
        end
    else
        if t < (initial+1)
            [val,idx] = max(T2(:,t));
            thetaall(:,:,t) = thetai{idx};
        elseif t >= (initial+1)
            [Tte(t),Itr(t)] = chartOC( T2(:,t),mT2,sd);
            Sall(:,:,t) = dall{Itr(t)};
            thetaall(:,:,t) = thetai{Itr(t)};
        end
    end
    timeforiteration = toc;
    formatSpec = 'Detect Time %d/%d with %3.3f seconds \n';
    fprintf(formatSpec,t,nT,timeforiteration)
end




if issave==1
    Snow = Sall;
    thetanow = thetaall;
else issave == 2
%     [~,Itr]= min(e);
    Snow = dall;
    thetanow = thetai;

end
    


end
