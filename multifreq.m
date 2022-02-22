%here code3a implemented in multifrequency
clear all;
addpath('C:\Local DiskD\Studies\Class Materials\L4T1\Thesis DOT\Paper\Exp5')
%addpath('C:\Local DiskD\Studies\Class Materials\L4T1\Thesis DOT\Paper\Exp3')

num=7;
pq=cell(1,num); % anom cascaded
Jac=cell(1,num); % J cascaded
data=cell(1,num);
ref=cell(1,num); %ref1 caccaded

freq=[ 300 500 600 400 700 200 100 ]; 

fwd_mesh=load_mesh('bckgN.rpt');
  pq{1}=load_data('anom2_data_300_n20.paa');
 pq{2}=load_data('anom2_data_500_n20.paa');
 pq{3}=load_data('anom2_data_600_n20.paa');
  pq{4}=load_data('anom2_data_400_n20.paa');
  pq{5}=load_data('anom2_data_700_n20.paa');
   pq{6}=load_data('anom2_data_200_n20.paa');
   pq{7}=load_data('anom2_data_100_n20.paa');
%   pq{8}=load_data('anom1_data_250_n20.paa');

recon_basis= [30 30];
iteration= 30;
lambda.value= 20;

% remove zeroed data
for i=1:num
pq{i}.paa(pq{i}.link(:,3)==0,:) = [];
data_link = pq{i}.link;

pq{i} = pq{i}.paa;
pq{i}(:,1) = log(pq{i}(:,1)); %take log of amplitude
pq{i}(:,2) = pq{i}(:,2)/180.0*pi; % phase is in radians and not degrees
pq{i}(pq{i}(:,2)<0,2) = pq{i}(pq{i}(:,2)<0,2) + (2*pi);
pq{i}(pq{i}(:,2)>(2*pi),2) = pq{i}(pq{i}(:,2)>(2*pi),2) - (2*pi);
pq{i} = reshape(pq{i}',size(pq{i},1)*2,1); 

end

fwd_mesh.link = data_link;
clear data


%*******************************************************
% Load or calculate second mesh for reconstruction basis

[fwd_mesh.fine2coarse,recon_mesh] = pixel_basis(recon_basis,fwd_mesh);

recon_mesh.type = fwd_mesh.type;
recon_mesh.link = fwd_mesh.link;
recon_mesh.source = fwd_mesh.source;
recon_mesh.meas = fwd_mesh.meas;
recon_mesh.dimension = fwd_mesh.dimension;

recon_mesh.element_area = ele_area_c(recon_mesh.nodes(:,1:2),...
                                         recon_mesh.elements);
                              
pixel.support = mesh_support(recon_mesh.nodes,...
                             recon_mesh.elements,...
                             recon_mesh.element_area);
                         
                         
% check for input regularization
    
    if size(pq{1},1)<2*size(recon_mesh.nodes,1)
        lambda.type = 'JJt';
    else
        lambda.type = 'JtJ';
    end
% Initiate projection error

pj_error = [];

for it = 1 : iteration
    
  it
  
  clear ref
  %**********
  % Calculate jacobian
  
  for i=1:num
      
  [Jac{i},data{i}]=jacobian_stnd(fwd_mesh,freq(i),recon_mesh);
  data{i}.amplitude(data_link(:,3)==0,:) = [];
  data{i}.phase(data_link(:,3)==0,:) = []; 
  
  Jac{i} = Jac{i}.complete;

  clear ref{i};
  
  ref{i}(:,1) = log(data{i}.amplitude);
  ref{i}(:,2) = data{i}.phase;
  ref{i}(:,2) = ref{i}(:,2)/180.0*pi;
  ref{i}(ref{i}(:,2)<0,2) = ref{i}(ref{i}(:,2)<0,2) + (2*pi);
  ref{i}(ref{i}(:,2)>(2*pi),2) = ref{i}(ref{i}(:,2)>(2*pi),2) - (2*pi);
  ref{i} = reshape(ref{i}',size(ref{i},1)*2,1);
   
  end
 
%   
  %*************
  clear anom ref1;
  anom=[];
  ref1=[];
  for i=1:num
    anom=[anom; pq{i}];
    ref1=[ref1; ref{i}];
  end
  data_diff = (anom-ref1);
  pj_error = [pj_error sum(abs(data_diff.^2))];
  pj_error(end)
  
   % Interpolate optical properties onto recon mesh
  [recon_mesh] = interpolatef2r(fwd_mesh,recon_mesh);
    
  % Normalize Jacobian wrt optical values
  N = [recon_mesh.kappa recon_mesh.mua];
  nn = length(recon_mesh.nodes);
  
  for j=1:num
  for i = 1 : nn
      Jac{j}(:,i) = Jac{j}(:,i).*N(i,1);
      Jac{j}(:,i+nn) = Jac{j}(:,i+nn).*N(i,2);
  end
  end
  
  clear nn N
  
  % build hessian
  clear J
  J=[];
  for i=1:num
  J=[J;Jac{i}];
  end
  
  [nrow,ncol]=size(J);
  if strcmp(lambda.type, 'JJt')
      Hess = zeros(nrow);
      Hess = (J*J');
      
      % regularising for amplitude and phase
      reg_amp = lambda.value*max(diag(Hess(1:2:end,1:2:end)));
      reg_phs = lambda.value*max(diag(Hess(2:2:end,2:2:end)));
      reg = ones(nrow,1);
      reg(1:2:end) = reg(1:2:end).*reg_amp;
      reg(2:2:end) = reg(2:2:end).*reg_phs;
      clear reg_*
  
      disp(['Amp Regularization        = ' num2str(reg(1,1))]);
      disp(['Phs Regularization        = ' num2str(reg(2,1))]);
      
      for i = 1 : nrow
          Hess(i,i) = Hess(i,i) + reg(i);
      end

      % Calculate update
      foo = J'*(Hess\data_diff);
      
  else %JTJ
 
      Hess = zeros(ncol);
      Hess = (J'*J);
      
      % regularising for amplitude and phase
      reg_kappa = lambda.value*max(diag(Hess(1:end/2,1:end/2)));
      reg_mua = lambda.value*max(diag(Hess(end/2+1:end,end/2+1:end)));
      reg = ones(ncol,1);
      reg(1:end/2) = reg(1:end/2).*reg_kappa;
      reg(end/2+1:end) = reg(end/2+1:end).*reg_mua;
      clear reg_*
      
       disp(['Kappa Regularization        = ' num2str(reg(1,1))]);
      disp(['Mua Regularization        = ' num2str(reg(end,1))]);

      % Add regularisation to diagonal - looped rather than creating a diaginal 
      % matrix as it is computational more efficient for large meshes
      for i = 1 : ncol
          Hess(i,i) = Hess(i,i) + reg(i);
      end

      % Calculate update
      foo = Hess\J'*data_diff;
  end
  
   % normalise back using optical parameters
  foo = foo.*[recon_mesh.kappa;recon_mesh.mua];
  
  % Update values
  recon_mesh.kappa = recon_mesh.kappa + (foo(1:end/2));
  recon_mesh.mua = recon_mesh.mua + (foo(end/2+1:end));
  recon_mesh.mus = (1./(3.*recon_mesh.kappa))-recon_mesh.mua;
  
  clear foo Hess Hess_norm tmp data_diff G
  
   %converge check
   if it ~= 1
    p = (pj_error(end-1)-pj_error(end))*100/pj_error(end-1);
    disp(['Projection error change   = ' num2str(p) '%']);
    if p <= 2 % stopping criteria is currently set at 2% decrease change
      disp('STOPPING CRITERIA REACHED');
     break
    end
  end


  % Interpolate optical properties to fine mesh
  [fwd_mesh,recon_mesh] = interpolatep2f(fwd_mesh,recon_mesh);
  
  % We dont like -ve mua or mus! so if this happens, terminate
  if (any(fwd_mesh.mua<0) | any(fwd_mesh.mus<0))
    disp('---------------------------------');
    disp('-ve mua or mus calculated...not saving solution');
  end
  
end

plotmesh(fwd_mesh);
%save_mesh(fwd_mesh,'C:\Local DiskD\Studies\Class Materials\L4T1\Thesis DOT\Paper\Exp5\anom2_7f_n20.rpt');
