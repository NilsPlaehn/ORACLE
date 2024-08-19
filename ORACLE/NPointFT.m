function Gp = NPointFT(MatInput,order,phi)
    N = numel(MatInput);
    tmp = 0;
    for j = 1:N
        tmp = MatInput(j).*exp(1i.*order.*phi(j))+tmp;
    end   
    Gp = 1./N.*tmp;
end