function varargout = montageArray(A,clim)

[m1,m2,m3] = size(A);

n = ceil(sqrt(m3));


M = zeros(m1*n, m2*n);

k=0;

for p=1:m1:(n*m1)
  for q=1:m2:(n*m2)
    k=k+1;
    if k>m3, 
      break
    end
    M(p:(p+m1-1),q:(q+m2-1)) = A(:,:,k);
  end
end

if nargout == 0
    if nargin == 1
        imagesc(M);
    else
        imagesc(M,clim);
    end   

    washold = ishold;
    hold on;
  
  [P,Q]= ndgrid(1:m1:((n+1)*m1),1:m2:((n+1)*m2));
  P = P -0.5;
  Q = Q -0.5;
  plot(Q,P,'color',get(gcf,'color'),'linewidth',3);
  plot(Q',P','color',get(gcf,'color'),'linewidth',3);
  plot(Q,P,'k','linewidth',1);
  plot(Q',P','k','linewidth',1);

  if ~washold,
      hold off;
  end;
  
end


if nargout == 1
  varargout{1} = M; 
end

  


