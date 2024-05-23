function out=RegulaCampo
global G
%G.C ->contem os out q tenho
%G.V ->contem os tam de cada out 
%keyboard
if G.tam>1
    G.camp=find(G.V==G.tam);
    if isempty(G.camp)
        G.control=0;
    else
        G.control=1;
        eval(['G.out = G.C' num2str(G.tam) ';'])
    end
    %    keyboard
    if G.control==0
        eval(['G.C' num2str(G.tam) '=campo2(G.location.x,G.location.y);'])
        eval(['G.out = G.C' num2str(G.tam) ';'])
        auu = G.V;
        G.V = zeros(1, length(auu)+1);
        G.V(1:end-1) = auu;
        G.V(end) = G.tam;
    end
    out.Br   = G.out.Br;
    out.Bz   = G.out.Bz;
    out.Bphi = G.out.Bphi;
    out.Ephi = G.out.Ephi;
    out.Bplx = G.out.Bplx;
    out.Bply = G.out.Bply;
    out.Leff = G.R0*G.B0./(0.001+(out.Br.^2 + out.Bz.^2).^0.5);
    out.a1= first_townsend_coeff(G.p0,out.Ephi,G.gas,0);
    %keyboard
else
    di=(0.5/2-G.location.x);
    s=real(exp(sqrt((di.^2+(G.location.y).^2))));
    if di==0
        out.Br=0;
        out.Bz=0;
    else
        alfa=atan(G.location.y/di);
        out.Br   = s*cos(alfa);
        out.Bz   = s*sin(alfa);
    end
    out.Bphi = G.R0*G.B0./(G.R0+G.location.x*G.a0);
    out.Ephi = -G.Vloop/2/pi./(G.R0+G.location.x*G.a0);
    out.Bplx = 0;
    out.Bply = 0;
    out.a1= first_townsend_coeff(G.p0,out.Ephi,G.gas,0);
    out.Leff = G.R0*G.B0./(0.001+(out.Br.^2 + out.Bz.^2).^0.5);
    %keyboard
end



end