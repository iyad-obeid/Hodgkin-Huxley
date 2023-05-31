%% Hodgkin-Huxley Fiber Implementation
% This is an implementation of a Hodgkin-Huxley fiber. The fiber is a
% concatenation of identical HH cells connected with a user-controlled
% conductance between cells
% It allows the user to provide a single stimulus at the time, strength,
% and duration of the user's choosing.

% Equations from Hodgkin, Huxley, J. Physiol (1952) 117, 500-544

% Iyad Obeid
% 25 November 2015
% Version 0.0

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
nCells = 100;

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
t = 0:dt:30;

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

tSTIM_START = 5;
tSTIM_DUR = 5;
STIM_STRENGTH = 20;

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
        if j == 1 % only stimulate the first cell
            if inRange( t(i) , tSTIM_START , tSTIM_START+tSTIM_DUR)
                Istim = STIM_STRENGTH;
            end
        end
        
        % DETERMINE CURRENT FROM CELL TO THE LEFT
        I_Left = 0;
        if ~isFirstCell(j)
            I_Left = gConnect*(v(j-1,i) - v(j,i));
        end
        
        % DETERMINE CURRENT FROM CELL TO THE RIGHT
        I_Right = 0;
        if ~isLastCell(j)
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

% define x-y cell locations
x = 1:100;
y = 0*x;

% define color map
clrs = hsv(130);

% define marker area
markerArea = 100;

% loop until ctrl-c
while(1)

    % define figure and make it pretty
    figure(2); clf;
    set(gcf,'color',[.2 .2 .2])

    % loop through time
    for i = 1:50:length(t)
        % draw the color-coded fiber
        scatter(x,y,markerArea,clrs(20+round(v(:,i)),:),'filled');
    
        % make the plot pretty
        set(gca,'visible','off');

        % update text
        text(50,0.5,sprintf('t=%.fs',t(i)),'color','w','fontsize',18);
        text(50,0.25,sprintf('CTRL-C to quit'),'color','w','fontsize',18);

        % force draw
        drawnow;
    end
    pause(1);
end


