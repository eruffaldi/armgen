
syms a alpha d q real
T = DH([a alpha, d q]);

syms oim1 doim1 oi doi qi dqi ddqi qim1 ddxi dxim1 ddxim1 dxm1  real
z0 = [0 0 1]';
oi = sym('oi_',[3,1]);
doi = sym('doi_',[3,1]);
DHP = T(1:3,4);
DHR = T(1:3,1:3);
oim1 = sym('oim1_',[3,1]);

% joint
Si = [];
Si.oi = DHR' *(oim1+dqi*z0);
Si.doi = DHR' *(doim1 + ddqi*z0 + dqi*cross(oim1,z0));
Si.dxi = DHR' * (dxim1 + cross(oi,DHP));
Si.ddxi = DHR' * (ddxim1 + cross(doi,DHP) + cross(oi,cross(oi, DHP)));

% Jacobian is recursive in the contribution so for every term we have a combination of the impacting qi

% joint no parent but moving
Si = [];
Sinpm.oi = DHR' *(dqi*z0);
Sinpm.doi = DHR' *(ddqi*z0);
Sinpm.dxi = DHR' * (cross(oi,DHP));
Sinpm.ddxi = DHR' * (cross(doi,DHP) + cross(oi,cross(oi, DHP)));

%
% D_qi(oi) = D_qi(DHR')(dqi*z0)
% D_dqi(oi) = DHR'*z0
% D_ddqi(oi) = 0
%
% D_qi(doi) = D_qi(DHR')*(ddqi*z0)
% D_dqi(doi) = 0
% D_ddqi(doi) = DHR'*z0
%
% then using the above ...
% D_qi(ddxi) = D_qi(DHR')*(...) + DHR'*(  D_qi(S(doi))*DHP + cross(doi,D_qi(DHP)) + D_qi(S(oi))S(oi)DHP + S(oi) D_qi(S(oi))*DHP + S(oi)S(oi)D_qi(DHP))
% D_dqi(ddxi) = DHR'*( D_dqi(S(oi))S(oi)DHP + S(oi) D_dqi(S(oi))*DHP)
% D_ddqi(ddxi) = DHR'*(  D_ddqi(S(doi))*DHP)
%
% When dealing with the higher level we have additional terms in the Jacobian to be summed AND
% also to take into account that D_dqi(doi) != 0 due to the added term, qhile D_ddqi(oi) is zero anyway
%
% J(oi) += DHR'(oim1)
% J(doi) += DHR'(doim1) + DHR'*dqi*cross(oim1,z0)
% J(ddx) += DHR' (ddxim1)
%
% The DHR' term contributes to all the D_qi
% In the doi we have a D_dqi term
% Finally each parent term contributes to the respectives higher level qim1
%
 


 