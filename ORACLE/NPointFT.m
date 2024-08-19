% Description: Discrete Fourier transform of bSSFP profile

% Author of function: 
% Nils MJ Plähn, Bern, Switzerland
% E-mail: nils.plaehn@students.unibe.ch
% Department of Diagnostic, Interventional and Pediatric Radiology (DIPR), Inselspital, Bern University Hospital, University of Bern, Switzerland
% Translation Imaging Center (TIC), Swiss Institute for Translational and Entrepreneurial Medicine, Bern, Switzerland



function Gp = NPointFT(MatInput,order,phi)
    N = numel(MatInput);
    tmp = 0;
    for j = 1:N
        tmp = MatInput(j).*exp(1i.*order.*phi(j))+tmp;
    end   
    Gp = 1./N.*tmp;
end
