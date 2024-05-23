function uinit = inicialfunc(location)
tam=length(location.x);%nrpontos;
uinit = zeros(Neq, tam); %condi inicial estados de u
uinit(1,:) = n0;              % u(1,:) -> n
uinit(4,:) = 1e-6;            % u(4,:) -> J_phi
uinit(2,:) = 0;               % u(2,:) -> J_x
uinit(3,:) = 0;               % u(3,:) -> J_y
uinit(5,:) = e*n0*Te0;         % u(5,:) -> pe
uinit(6,:) = e*n0*Te0;         % u(6,:) -> pion
uinit(7,:) = 7.89e11*p0;      % u(7,:) -> v_en
uinit(8,:) = 0;               % u(8,:) -> v_in
uinit(9,:) = out1.v_ei;       % u(9,:) -> v_ei
uinit(10,:) = 0;              % u(10,:) -> v_loss
setInitialConditions(model,uinit); 
% for j=1:1:Neq %checando as derivadas em x das Neq equações
%     if state.time==0
%         state.ux(j,:)=0; state.uy(j,:)=0;
%     end
%     if isnan(state.ux(j,1))
%         state.ux(j,:)=0;
%     end
%     if isnan(state.uy(j,1))
%         state.uy(j,:)=0;
%     end
% end

% if isnan(state.u(1,1))
%     state.u(1,:)=n0;
% end
% if isnan(state.u(3,1))
%     auxPe=e*Te0*n0*ones(1,tam);
% else
%     auxPe=state.u(3,:);
% end
% 
% if isnan(state.u(4,1))
%     auxPi=e*Te0*n0*ones(1,tam);
% else
%     auxPi=state.u(4,:);
% end
% 
% for j=1:1:Neq %checando as derivadas em x das Neq equações
%     
%     if isnan(state.ux(j,1))
%         state.ux(j,:)=0;
%     end
%     if isnan(state.uy(j,1))
%         state.uy(j,:)=0;
%     end
% end