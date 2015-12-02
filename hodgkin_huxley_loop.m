%% Hodgkin-Huxley Loop Implementation
% This is an implementation of a Hodgkin-Huxley loop. The loop is a
% concatenation of identical HH cells connected with a user-controlled
% conductance between cells.
% It allows the user to provide two stimuli at the time, strength,
% and duration of the user's choosing.
% The default stimuli establish an 'infinite' loop of activation

% Equations from Hodgkin, Huxley, J. Physiol (1952) 117, 500-544

% Iyad Obeid
% 25 November 2015
% Version 0.1

clear;
figure(1);clf;

%%%%%%%%%%%%%%%%%%
% DEFINE CONSTANTS
%%%%%%%%%%%%%%%%%%

% Nernst Potentials
Ena = 115; Ek = -12; El = 10.6;

% Maximum Conductances
gna = 120; gk = 36; gl = 0.3;

% Membrane Capacitance
C = 1;

% Number of Cells
nCells = 500;

% Connectivity between Cells
gConnect = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE VOLTAGE-DEPENDENT GATE ACTIVATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Likelihoods of gates opening
an = @(u) (0.1-0.01*u)/(exp(1-0.1*u)-1);
am = @(u) (2.5-0.1*u)/(exp(2.5-0.1*u)-1);
ah = @(u) 0.07*exp(-u/20);

% Likelihoods of gates closing
bn = @(u) 0.125*exp(-u/80);
bm = @(u) 4*exp(-u/18);
bh = @(u) 1/(exp(3-0.1*u)+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE FORMULAE FOR STEADY STATE GATE ACTIVATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m_inf = @(u) am(u) / ( am(u) + bm(u) );
n_inf = @(u) an(u) / ( an(u) + bn(u) );
h_inf = @(u) ah(u) / ( ah(u) + bh(u) );

%%%%%%%%%%%%%
% DEFINE TIME
%%%%%%%%%%%%%

dt = 0.01;
t = 0:dt:100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE STATE VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = zeros(nCells,length(t)) + m_inf(0);
n = zeros(nCells,length(t)) + n_inf(0);
h = zeros(nCells,length(t)) + h_inf(0);
v = zeros(nCells,length(t));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE STIMULUS STRENGTH, DURATION, & DELAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stims(1).tSTIM_START    = 0;
stims(1).tSTIM_DUR      = 5;
stims(1).STIM_STRENGTH  = 50;
stims(1).cell           = 1;

stims(2).tSTIM_START    = 33;
stims(2).tSTIM_DUR      = 9;
stims(2).STIM_STRENGTH  = 5000;
stims(2).cell           = round(nCells/4);

%%%%%%%%%%%%%
% DEFINE MISC
%%%%%%%%%%%%%

inRange = @(x,a,b) (x>=a) & (x<b);
isFirstCell = @(j) j == 1;
isLastCell  = @(j) j == nCells;

%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%

for i = 1:length(t)-1
    if mod(i,100)==0
        fprintf('%i\n',i);
    end
    
    for j = 1:nCells
        
        % EXTRACT MEMBRANE VOLTAGE
        u = v(j,i);
        
        % SOLVE FOR MEMBRANE CURRENTS
        Ik  = gk  * n(j,i)^4           * ( u - Ek  );
        Ina = gna * m(j,i)^3 * h(j,i)  * ( u - Ena );
        Il  = gl  *                      ( u - El  );
        I_mem = Ik + Ina + Il;
        
        % DETERMINE STIMULUS CURRENT, IF ANY
        Istim = 0;
        for k = 1:length(stims)
            if (stims(k).cell==j) && inRange( t(i) , stims(k).tSTIM_START , stims(k).tSTIM_START+stims(k).tSTIM_DUR)
                Istim = stims(k).STIM_STRENGTH;
            end
        end
        
        % DETERMINE CURRENT FROM CELL TO THE LEFT
        if isFirstCell(j)
            I_Left = gConnect*(v(nCells,i) - v(j,i));
        else
            I_Left = gConnect*(v(j-1,i) - v(j,i));
        end
        
        % DETERMINE CURRENT FROM CELL TO THE RIGHT
        if isLastCell(j)
            I_Right = gConnect*(v(1,i) - v(j,i));
        else
            I_Right = gConnect*(v(j+1,i) - v(j,i));
        end
        
        
        % DEFINE THE STATE VARIABLE DERIVATIVES
        dv = (Istim + I_Left + I_Right - I_mem)/C;
        dm = am(u) * (1-m(j,i)) - bm(u) * m(j,i);
        dh = ah(u) * (1-h(j,i)) - bh(u) * h(j,i);
        dn = an(u) * (1-n(j,i)) - bn(u) * n(j,i);
        
        % USE FORWARD EULER TO INCREMENT THE STATE VARIABLES
        v(j,i+1) = v(j,i) + dv*dt;
        m(j,i+1) = m(j,i) + dm*dt;
        h(j,i+1) = h(j,i) + dh*dt;
        n(j,i+1) = n(j,i) + dn*dt;
        
    end
end


% plot the results
figure(1);clf;
plot(t,v( 1:10:end ,:)); % plot the voltage at every 10th cell
axis([t(1) t(end) -20 120]);
xlabel('time (ms)');
ylabel('membrane voltage (mV)');

%%%%%%%%%%%%%%%%%%%%%
% ANIMATE THE RESULTS
%%%%%%%%%%%%%%%%%%%%%

% define the colormap
clrs = hsv(130);

% define x-y coordinates for cell locations
x = cos(2*pi*(1:nCells)/nCells-pi/2);
y = sin(2*pi*(1:nCells)/nCells+pi/2);

markerArea = 100;

% just keep on looping!
while(1) % hit control-c to quit
    
    % initialize and color the figure
    figure(2); clf;
    set(gcf,'color',[.2 .2 .2])
    
    % animate every 100th time-step
    for i = 1:100:length(t)
        
        % draw the color-coordinated cells
        % the '20' and '130' are custom offsets to get the colors looking
        % right
        scatter(x,y, ...
            markerArea, ...
            clrs( min(20+round(v(:,i)),130), :), ...
            'filled'...
            );
        
        % add the time text
        text(-0.2,0,sprintf('t=%.1fs',t(i)),'color','w','fontsize',18);
        text(-0.2,-0.2,sprintf('CTRL-C to quit'),'color','w','fontsize',18);
        % make it pretty
        axis equal
        set(gca,'visible','off');
        
        % force re-draw
        drawnow;
        
    end
    
    pause(1);
end
