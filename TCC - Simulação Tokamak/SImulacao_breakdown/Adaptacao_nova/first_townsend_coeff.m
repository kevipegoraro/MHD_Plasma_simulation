function a = first_townsend_coeff(p,E,gas,iplot)

switch gas
    case 'H2'
        A = 4.8; % In units of [ 1 / ( cm . Torr ) ]
        B = 136; % In units of [ V / ( cm . Torr ) ]
    case 'N2'
        A = 11.8;
        B = 325;
    case 'O2'
        A = 6.5;
        B = 190;
    case 'He'
        A = 2.8;
        B = 77;
    case 'Ne'
        A = 4.4;
        B = 111;
    case 'Ar'
        A = 11.5;
        B = 176;
    case 'Kr'
        A = 15.6;
        B = 220;
    case 'Xe'
        A = 24;
        B = 330;
    case 'CH4'
        A = 17;
        B = 300;
    case 'CF4'
        A = 11;
        B = 213;
end

% Converting units from A[ 1 / ( cm . Torr ) ] into A[ 1 / ( m . Pa) ]
% and from B[ V / ( cm . Torr ) ] into A[ V / ( m . Pa) ]
A = A/1e-2/133.322;
B = B/1e-2/133.322;

% First Townsend coefficient
a = A*p.*exp(-B*p./E);

% Plotting
if iplot
    figure
    plot(E,a)
    xlabel('Electric Field ( V / m )')
    title(['First Townsend Coefficient ( 1 / m ) for ' gas])

    figure
    plot(E,p)
    xlabel('Gas Pressure ( Pa )')
    title(['First Townsend Coefficient ( 1 / m ) for ' gas])
end

end      