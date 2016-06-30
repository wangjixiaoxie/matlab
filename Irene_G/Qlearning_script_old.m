% Actions cost energy and are thus associated with a -0.2 reward.

% Alternating correctly between the two outside arms results in a reward of

% 1-0.2. Impossible actions have reward 0.

 

% e is exploration factor

e= 0;

 

% Transition probability array (SxSxA)

P(:,:,1)=[0,1-e,0,0;0,0,1-e,0;0,0,0,0;0,0,1-e,0];

P(:,:,2)=[0,0,0,0;1-e,0,0,0;0,0,0,1-e;1-e,0,0,0];

P(:,:,3)=[0,0,1-e,0;0,0,0,0;0,0,0,0;0,0,0,0];

P(:,:,4)=[0,0,0,0;0,0,0,0;1-e,0,0,0;0,0,0,0];

P(:,:,5)=[1-e,0,0,0;0,1-e,0,0;0,0,1-e,0;0,0,0,1-e];

 

% Constructing reward matrix (SxSXA)

 

go1Right = [0,0.8,0,0;0,0,0.8,0;...

             0,0,0,0;0,0,-0.2,0];

         

go1Left = [0,0,0,0;-0.2,0,0,0;...

            0,0,0,0.8;0.8,0,0,0];

        

go2Right = [0,0,-0.2,0;0,0,0,0;...

            0,0,0,0;0,0,0,0];

        

go2Left = [0,0,0,0;0,0,0,0;...

            -0.2,0,0,0;0,0,0,0];

 

goSame =   [-0.2,0,0,0;0,-0.2,0,0;...

            0,0,-0.2,0;0,0,0,-0.2];

   

        

R(:,:,1)=(go1Right);

R(:,:,2)=(go1Left);

R(:,:,3)=(go2Right);

R(:,:,4)=(go2Left);

R(:,:,5)=(goSame);

 

% choose temperature, value btw 0 and 1, 0= no noise

t=0.4;

 

% discount factor, represents how far into the future algorith looks, 0-1

discount=0.9;

 

% choose trial number

Nr=300;

 

% mdp_Q_learning   Evaluation of the matrix Q, using the Q learning algorithm 

%

% Arguments

% -------------------------------------------------------------------------

% Let S = number of states, A = number of actions

%   P(SxSxA)  = transition matrix 

%              P could be an array with 3 dimensions or 

%              a cell array (1xA), each cell containing a sparse matrix (SxS)

%   R(SxSxA) or (SxA) = reward matrix

%              R could be an array with 3 dimensions (SxSxA) or 

%              a cell array (1xA), each cell containing a sparse matrix (SxS) or

%              a 2D array(SxA) possibly sparse  

%   discount  = discount rate in ]0; 1[

%   N(optional) = number of iterations to execute, default value: 10000.

%                 It is an integer greater than the default value. 

% Evaluation --------------------------------------------------------------

%   Q(SxA) = learned Q matrix 

%   V(S)   = learned value function.

%   policy(S) = learned optimal policy.

%   mean_discrepancy(N/100) = vector of V discrepancy mean over 100 iterations

%             Then the length of this vector for the default value of N is 100.

 

 

 

 

% check of arguments

if (discount <= 0 || discount >= 1)

    disp('--------------------------------------------------------')

    disp('MDP Toolbox ERROR: Discount rate must be in ]0,1[')

    disp('--------------------------------------------------------')   

elseif (nargin >= 4) && (N < 10000)

    disp('--------------------------------------------------------')

    disp('MDP Toolbox ERROR: N must be upper than 10000')

    disp('--------------------------------------------------------') 

else

 

    % initialization of optional arguments

    % I set the trial number here? Originally it was 10000

N    if (nargin < 4); N=Nr; 
end;      

    

    % Find number of states and actions

    if iscell(P)

        S = size(P{1},1);

        A = length(P);

    else

        S = size(P,1);

        A = size(P,3); 

    end;

    

    % Initialisations

    Q = zeros(S,A);

    dQ = zeros(S,A);

    mean_discrepancy = [];

    discrepancy = [];

 

    % Initial state choice

    s = randi([1,S]);

    

    % Initialize arrays for probabilities of right action for each state

    ProbAll=[];

    Prob01All=[]; % when left, going one right=  inbound

    Prob22All=[]; % when right, going one left= inbound

    Prob12All=[]; % when middle after left going one right outbound

    Prob32All=[]; % when middle after right going one left outbound

 

    for n=1:N

 

        % Reinitialisation of trajectories every 100 transitions

        % Animal runs about 100 trials a day, so should be alright

        if (mod(n,100)==0); s = randi([1,S]); end;

        

        % Action choice : greedy with increasing probability

        % probability 1-(1/log(n+2)) can be changed, this is some e factor, which

        % should be ok

        pn = rand(1);

        if (pn < (1-(1/log(n+2))))

          [nil,a] = max(Q(s,:));

        else

          a = randi([1,A]);

        end;

 

        % Simulating next state s_new and reward associated to <s,s_new,a> 

        p_s_new = rand(1);

        p = 0; 

        s_new = 0;

        while ((p < p_s_new) && (s_new < S)) 

            s_new = s_new+1;

            if iscell(P)

                p = p + P{a}(s,s_new);

            else   

                p = p + P(s,s_new,a);

            end;

        end; 

        if iscell(R)

            r = R{a}(s,s_new); 

        elseif ndims(R) == 3

            r = R(s,s_new,a); 

        else

            r = R(s,a); 

        end;

 

        % Updating the value of Q   

        % Decaying update coefficient (1/sqrt(n+2)) can be changed

     

        delta = r + discount*max(Q(s_new,:)) - Q(s,a);

        dQ = (1/sqrt(n+2))*delta;

        Q(s,a) = Q(s,a) + dQ;

       

        

        % Current state is updated

        s = s_new;

        

        % Calculate probability of Q(s,a)using Boltzmann distribution

        

Prob1run=[];

 

for ii=1:S

    

Prob=(exp(Q(ii,:))/t)/sum(exp(Q(ii,:)/t));

Prob1run=[Prob1run; Prob];

 

end

 

% Select prob of right action for each state in 1 run

Prob01=Prob1run(1,1); 

Prob22=Prob1run(3,2);

Prob12=Prob1run(2,2);

Prob32=Prob1run(4,2);

 

        % Computing and saving maximal values of the Q variation  

        discrepancy(mod(n,100)+1) = abs(dQ);  

    

        % Computing means all over maximal Q variations values  

        if (length(discrepancy) == 100)     

           mean_discrepancy = [ mean_discrepancy mean(discrepancy)];

           discrepancy = [];

        end;  

        

% Save probabilites for right action in each state for each run 

ProbAll=[ProbAll;Prob1run];

Prob01All=[Prob01All;Prob1run(1,1)];

Prob22All=[Prob22All;Prob1run(3,2)]; 

Prob12All=[Prob12All;Prob1run(2,2)]; 

Prob32All=[Prob32All;Prob1run(4,2)];

    end;

 

    % Compute the value function and the policy

    [V, policy] = max(Q,[],2);        

 

end;

 

% Plot figure: green inbound trials, blue outbound trials

 

figure;

x=1:1:Nr;

plot(x,Prob01All,'g',x,Prob22All,'g--',x,Prob12All,'b',x,Prob32All,'b--');

xlabel ('Trial Number');

ylabel ('Probability of right action given current state');