function [VI, NMI] = VI(U,V)

n = length(U);

Ku = length(unique(U));
Kv = length(unique(V));

Du = zeros(Ku,1);
for i = 1:Ku
    Du(i) = sum(U== i);
end
Du = (1/n).*Du;

Dv = zeros(Kv,1);
for i = 1:Kv
    Dv(i) = sum(V == i);
end
Dv = (1/n).*Dv;


JD = Du*(Dv');

Hu = -sum(Du.*log(Du));
Hv = -sum(Dv.*log(Dv));
I = nansum(log(JD./(Du*Dv')).*JD,'all');

VI = Hv + Hu-2*I;

NMI = sqrt(I^2/(Hv*Hu));