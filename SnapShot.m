tic
                      %%% Korosh Mahmoodi 021021
clc;  
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% for Bottom-up version set: 
            %%% TimeShotDown = Trial + 1 (means top-down govern will not happen)
            %%% TimeOuTiSolaTion = Trial + 1
            %%% InitialInfected =  1 
            %%% NumberIsolated = 0
            
            
            %%%%%%%%% For Top-Down version set:
            %%% OthersSympt = 0
            %%% SelfSympt = 0 
            %%% IsolationPeriod = TimeOuTiSolaTion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                            Trial = 2e4 ;  %%% number of itteration
                            
                                   TimeOutforPairSearch = Trial ; % 1e4 ; %% limit for pair search
                                                                                     
                                     Wealth = zeros( Trial , 1) ;                                              
                                     Delay = zeros( Trial , 1) ;
                                     WealthReal = zeros( Trial , 1) ; 
                                                
                                     Meet12  = zeros(2*Trial,1);

 mov = VideoWriter('SA.avi') ;
open(mov)
                                     
AccumPairWealth = zeros(Trial,1);                                                
                            
                                                             ENSRatio = 1 ;
                                                                                                                    
%                                                                       ImitIsol = 1 ;
                                                          Size = 10 ; %%% number of agents
                       AccumAllWealth = zeros(Trial,Size ) ;                                                
        
                            InitialInfected =  4 ;  %%% number of initial infecteds
                            NumberIsolated = 0 ; %% round(  YYYY(ss13) *  Size  ) ;  %% for top down version
                            
                            TimeVirus = 0.5 * Trial ; %%% time when one virus infects one agent (or a ratio of agents in the top-down version)
                            LifeTimeVirus = 2e5 ; %1e4 ;  

                            TimeShotDown = Trial +1 ; %%% TimeVirus + 0.5 * LifeTimeVirus ; %%%  TimeVirus   +  XXXX(rr13) * LifeTimeVirus  ;  %%%
                            TimeOuTiSolaTion = Trial +1 ; %%% TimeShotDown + 0.5 * LifeTimeVirus ; %%% TimeShotDown   +  2 * LifeTimeVirus  ;  %%%

                                                        SymptomPeriod = 100 ; % 0.5 * LifeTimeVirus  ;  %%%% Taking tests lowers this
                                                            IsolationPeriod =  2 * LifeTimeVirus ;  %%% fixed
                                          %%%%% more isolation is useless bec virus stikes 
                                          %%%%%  less than 0.5 means agents cant isolate more need go out food
                                          %%%% or a distribution of isolation times

                                                        DeadChance  = 0.1 ; % 0.2 ;  %% chance that infected dies after LifeTime of virus
                                                        %%% ventilator etc.
                                                        InfectionChance = 0.2 ; % 0.3 ; %% chance that infection spreads to the partner
                                                         %%% mask and
                                                         %%% washing
                                                         
                                                        OthersSympt = 0 ; % 0.25 ;%%% XXXX(rr13) ;  %%% chance that an agent goes to isolation seeing symptom in its pair
                                                        SelfSympt = 0 ; %0.25 ;%%%0.25 ;%%%XXXX(rr13) ; %% chance that agent goes to isolation if sees infection in itself

                                     AveIm = zeros(Trial, 1) ;

                                      MeanDied = zeros( Trial, 1 ) ;
                                      MeannInfected  = zeros( Trial, 1 ) ;
                                      MeanVaccinated = zeros( Trial, 1 ) ;
                                          MeanCC = zeros( Trial, 1 ) ;
                                          MeanCDDC = zeros( Trial, 1 ) ;
                                          MeanDD = zeros( Trial, 1 ) ;

                             MeannIsolatedInfec= zeros( Trial, 1 ) ;
                             MeannIsolatedHealthy= zeros( Trial, 1 ) ;
                             MeannFREEInfect= zeros( Trial, 1 ) ;
                             MeannFREEHealthy= zeros( Trial, 1 ) ;
                                         Meanfield = zeros( Trial, 1 ) ;
                             
for zz = 1 : ENSRatio
zz
                                             Died = zeros( Size , 1) ;
                                             Infected = zeros( Size , 1) ;
                                             Vaccinated = zeros( Size , 1) ;
                                             Isolated = zeros( Size , 1) ;
                                             
                                             TimeInfected = zeros( Size , 1) ;
                                             TimeIsolated = zeros( Size , 1) ;
                                             
                                        %%%%%%%%%%%% SA

                           % In SA, there are three decisions for each agent to make:
% not-Connect or Connect (1 or 2)
% Defect or Cooperate (1 or 2)
% not-Trust or Trsut (1 or 2)
% For each choice there is a corresponding adaptive propensity/probability

% Constants relating payoff (as feedback) to the change of the propensity/probability of the corresponding decision made
Chi_Connection = 0.1 ; % SAC mechanism (connect decision)
Chi_RL = 0.1 ; %  SAL mechanism (Cooperation or Defection decision)
Chi_Trust = 0.1 ; % SAT mechanism (trust decision)
% You can deactivate SAT, SAC, or both mechanisms by seting the corresponding Chi constant(s) to zero (for deactivating SAT you must set trust initial decision and propensity to 1)

% Parameters of the Prisoner's Dilema game's
s = 0 ; % s is the payoff of the cooperator agent if another agent defected
p = 0 ; % p is the payoff of the agents if both defected
tc = 0.9 ; % temptation to cheat. (1+ tc) is the payoff of the defector agent if the other agent cooperated (tc < s + 1 and  tc < 1)
                
DemoTime = 1*Trial ; % Time for the demo

%%%%%%%%%%%%%%%%%%%% initial conditions
Ratio_CC = zeros(Trial, 1) ; CC = 0 ; % Ratio of Cooperation-Cooperation decisions
Ratio_TNT = zeros(Trial, 1) ; TNT = 0 ; % Ratio of Trust-Not Trust decisions

Out1 = zeros(Size, 1) ; % Previous payoff of agents
Out2 = zeros(Size, 1) ; % Current payoff of agents

Connection_decision = zeros(Size, 1) ; % Not-connct as 1 and connect as 2
D_C_decision = zeros(Size, 1) ; % D (defection) as 1 and C (cooperation) as 2
Trust_decision = zeros(Size, 1) ; % Not-Trust as 1 and Trust as 2

D_C_decisionColor = zeros(Trial, Size) ; % Assigns color for the nodes at each trial

Infected_decisionColor = zeros(Trial, Size) ;
Isolated_decisionColor = zeros(Trial, Size) ;
Vaccinated_decisionColor = zeros(Trial, Size) ;
Died_decisionColor = zeros(Trial, Size) ;

P_Connection =  zeros(Size, Size -1) ; % Tendensy of the agents to connect/play
P_RL = zeros(Size, Size -1, 2) ; % P_RL and 1-P_RL are the propensity of decison D and C, respectively
P_Trust = zeros(Size,Size -1, 2) ; % P_Trust and 1-P_Trust are the propensity of decison "not to trust" and "trust" the decison of another agent, respectively

Conect2ty = zeros(Trial, Size, Size) ; % Records the Connection matrix of each trial


% Initial conditions
for jjj = 1 : Size
    r = rand ;
    if r < 0.5
        D_C_decision(jjj) = 1 ;
    else
        D_C_decision(jjj) = 2 ;
    end

    r = rand ;
    if r < 0.5
        Trust_decision(jjj) = 1 ;
    else
        Trust_decision(jjj) = 2 ; % = 1 if you want to eliminate the SAT mechinisem. Also in line 78 use P_Trust(i , j,  1) =  1 ;
    end
end

for i = 1 : Size
    for j = 1 : Size
        for k = 1 : 2 % number of choices/propensities
            if i ~= j
                P_RL(i , j, k) = 0.5 ;
                P_Trust(i , j,  k) =  0.5 ;
            end
        end
    end
end

for i = 1 : Size
    for j = 1 : Size - 1 % number of choices for each agent to connect
        P_Connection(i, j) = 1/(Size - 1) ;
    end
end


AveTrust = zeros(Trial, 1) ;

C1C1=0;

gh =  1 ;


for ti =  2 : Trial
ti/Trial
    P_Connection_Used = ones(Size, 1) ; % Tracks if the Connection decision mechanism has used (2) in the final decision making or not (1)
    P_RL_Used = ones(Size, 1) ; % Tracks if the RL decision mechanism has used (2) in the final decision making or not (1)
    P_Trust_Used = ones(Size, 1) ; % Tracks if the Trust decision mechanism has used (2) in the final decision making or not (1)
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIAL INFECTED

if ti == TimeVirus  %%% in top down version a ratio randomly gets infected

L0 = InitialInfected ;
out = randperm(Size,L0)  ;

for oo91 = 1 : InitialInfected
       Infected(out(oo91)) = 1 ;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ISOLATED AGENTS

if ti == TimeShotDown  %%% in top down a ration randomly gets isolated
    
n0 = Size  ;
L0 = NumberIsolated ;
out0 = randperm(n0,L0)  ;

for oo11 = 1 : NumberIsolated
       Isolated(out0(oo11)) = 1 ;
end
end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Picking agents to play    

 Surpass = 0 ;
      m = floor( 1 + Size * rand) ;
     tiDelay = 1 ;        
       while Isolated(m) == 1  || Died(m) == 1                 
      
          m = floor( 1 + Size * rand) ;
          tiDelay =  tiDelay + 1 ;    

if tiDelay >= TimeOutforPairSearch   && ( Died(m) == 0   )  
    
    Isolated(m) = 0 ;
    n = floor( 1 + Size * rand) ;
    tiDelay =  tiDelay + 1 ;    

    while m == n  ||  Isolated(n) == 1   || Died(n) == 1         
   
          n = floor( 1 + Size * rand) ;
          tiDelay =  tiDelay + 1 ;    

if  Died(n) == 0   &&   n~=m       
    
    Isolated(n) = 0 ;
    Surpass = 1 ;
end     
   end
end
      end
      
 if Surpass == 0     
 
         n = floor( 1 + Size * rand) ;
         tiDelay = tiDelay + 1 ;  
      while m == n                   ||  Isolated(n) == 1  || Died(n) == 1         
          n = floor( 1 + Size * rand) ;
             tiDelay =  tiDelay + 1 ;   

if tiDelay > TimeOutforPairSearch   && (  Died(n) == 0 )   &&   n~=m   
    
    Isolated(n) = 0 ;
    Surpass = 1 ;    
end     
      end
 end     
 
            
while Connection_decision(n) ~= m || Connection_decision(m) ~= n
            
                                      %%%%%%%% should forced to out of isolation be back 

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Picking units to play    
      m = floor( 1 + Size * rand) ;
      tiDelay = tiDelay + 1 ; 

       while Isolated(m) == 1  || Died(m) == 1                 
      
          m = floor( 1 + Size * rand) ;
          tiDelay =  tiDelay + 1 ;    

if tiDelay >= TimeOutforPairSearch   &&  Died(m) == 0   
    
    Isolated(m) = 0 ;
   
    n = floor( 1 + Size * rand) ;
    tiDelay =  tiDelay + 1 ;    

    while m == n   ||   Isolated(n) == 1   || Died(n) == 1         
   
          n = floor( 1 + Size * rand) ;
          tiDelay =  tiDelay + 1 ;    

if    Died(n) == 0   &&   n~=m       
    Isolated(n) = 0 ;
    Surpass = 1 ;
end     
   end
end
      end
      
      
 if Surpass == 0     
 
         n = floor( 1 + Size * rand) ;
         tiDelay = tiDelay + 1 ;  
      while m == n                   || Isolated(n) == 1 || Died(n) == 1         
          n = floor( 1 + Size * rand) ;
          tiDelay =  tiDelay + 1 ;   

if tiDelay > TimeOutforPairSearch   &&  Died(n) == 0   &&     n~=m    
    
    Isolated(n) = 0 ;
    Surpass = 1 ;    
end     
      end
 end     



        PCum_Connection = cumsum(P_Connection(m , :) ) ;
        r = rand ;
        for i = 1 : Size - 1
            if  i < m
                if r < PCum_Connection(i)
                    Connection_decision(m) = i ;
                    break
                end
            end
            if  i > m
                if r < PCum_Connection(i)
                    Connection_decision(m) = i + 1 ;
                    break
                end
            end
        end

        PCun_Connection = cumsum(P_Connection(n , :) ) ;
        r = rand ;
        for i = 1 : Size - 1
            if  i < n
                if r < PCun_Connection(1,i)
                    Connection_decision(n) = i ;
                    break
                end
            end
            if  i > n
                if r < PCun_Connection(1,i)
                    Connection_decision(n) = i + 1 ;
                    break
                end
            end
        end
 
end 


P_Connection_Used(m) = 2 ; % Connection decision mechanism used
    P_Connection_Used(n) = 2 ;

%%%%%%%%%%
                          Delay(ti) = tiDelay  ;  %%% number of attempts to pair

    % Decision of C or D
    % D_C_decision(m) = Decision(P_RL(m, n, :)) ;
    % D_C_decision(n) = Decision(P_RL(n, m, :)) ;

PCum_RL = cumsum(P_RL(m, n, :) ) ;
    r = rand ;
    for i = 1 : 2
        if r < PCum_RL(1,1,i)
            D_C_decision(m) = i ;
            break
        end
    end

    PCun_RL = cumsum(P_RL(n, m, :) ) ;
    r = rand ;
    for i = 1 : 2
        if r < PCun_RL(1,1,i)
            D_C_decision(n) = i ;
            break
        end
    end

    P_RL_Used(m) = 2 ; % i.e. RL decision mechanism has used
    P_RL_Used(n) = 2 ;

    %  Decision to not-trust or trust the decision of the other agent
    % Trust_decision(m) = Decision(P_Trust(m, n, :)) ;
    % Trust_decision(n) = Decision(P_Trust(n, m, :)) ;
    PCum_Trust = cumsum(P_Trust(m , n, :) ) ;
    r = rand ;
    for i = 1 : 2
        if r < PCum_Trust(1, 1, i)
            Trust_decision(m) = i ;
            break
        end
    end

    PCun_Trust = cumsum(P_Trust(n ,m , :) ) ;
    r = rand ;
    for i = 1 : 2
        if r < PCun_Trust(1,1,i)
            Trust_decision(n) = i ;
            break
        end
    end

    D_C_m = D_C_decision(m) ;
    D_C_n = D_C_decision(n) ;

    if  Trust_decision(m) == 2
        D_C_decision(m) = D_C_n ;

        P_Trust_Used(m) = 2 ; % Trust decision mechanism used
        P_RL_Used(m) = 1 ; % RL decision mechanism eliminated
    end
    if  Trust_decision(n) == 2
        D_C_decision(n) = D_C_m ;

        P_Trust_Used(n) = 2 ;
        P_RL_Used(n) = 1 ;
    end

    if  Trust_decision(m) * Trust_decision(n) == 2
        TNT = TNT + 1 ;
    end
    Ratio_TNT(ti) = Ratio_TNT(ti) + TNT/ ti ;


    % Payoffs from the Prisoner's Dilemma game
    % Out2(m) = Payoff_PD(D_C_decision(m), D_C_decision(n), s, p, tc) ;
    % Out2(n) = Payoff_PD(D_C_decision(n), D_C_decision(m), s, p, tc) ;
    if D_C_decision(m) == 2
        ggn = 1 ;
    else
        ggn = 0 ;
    end

    if D_C_decision(n) == 2
        ggm = 1 ;
    else
        ggm = 0 ;
    end

    if D_C_decision(m) == 2
        Out2(m) = ggm * (1) + (1 - ggm) * (-s) ;
    else
        Out2(m) = ggm * (1 + tc) + (1 - ggm) * (p) ;
    end

    if D_C_decision(n) == 2
        Out2(n) = ggn * (1) + (1 - ggn) * (-s) ;
    else
        Out2(n) = ggn * (1 + tc) + (1 - ggn) * (p) ;
    end


    % Updating the propensities
    PConnectionm = zeros(1, Size-1);
    PConnectionn = zeros(1, Size-1);
    for k = 1 : Size-1
        PConnectionm(k) = P_Connection(m, k) ;
        PConnectionn(k) = P_Connection(n, k) ;
    end

    PRLm = zeros(1, 2) ;
    PRLn = zeros(1, 2) ;
    PTrustm = zeros(1, 2) ;
    PTrustn = zeros(1, 2) ;
    for k = 1 : 2

        PRLm(k) = P_RL(m, n, k) ;
        PRLn(k) = P_RL(n, m, k) ;

        PTrustm(k) = P_Trust(m, n, k) ;
        PTrustn(k) = P_Trust(n, m, k) ;
    end

     if Connection_decision(m) < m
        Cdecisionm = Connection_decision(m) ;
    else
        Cdecisionm = Connection_decision(m) - 1 ;
    end
    if Connection_decision(n) < n
        Cdecisionn = Connection_decision(n) ;
    else
        Cdecisionn = Connection_decision(n) - 1 ;
    end

    P_Connection(m, :) = SA_Update(P_Connection_Used(m), PConnectionm, Cdecisionm, Out2(m) , Out1(m) , Chi_Connection) ;
    P_Connection(n, :) = SA_Update(P_Connection_Used(n), PConnectionn, Cdecisionn, Out2(n) , Out1(n) , Chi_Connection) ;

    P_RL(m, n, :) = SA_Update(P_RL_Used(m), PRLm , D_C_decision(m), Out2(m), Out1(m), Chi_RL) ;
    P_RL(n, m, :) = SA_Update(P_RL_Used(n), PRLn , D_C_decision(n), Out2(n), Out1(n), Chi_RL) ;

    P_Trust(m, n, :) = SA_Update(P_Trust_Used(m), PTrustm, Trust_decision(m) , Out2(m), Out1(m), Chi_Trust) ;
    P_Trust(n, m, :) = SA_Update(P_Trust_Used(n), PTrustn, Trust_decision(n) , Out2(n), Out1(n), Chi_Trust) ;

    % Updating the previous payoffs
    Out1(m) = Out2(m) ;
    Out1(n) = Out2(n) ;

    if  D_C_decision(m) == 2  &&  D_C_decision(n) == 2
        CC = CC + 1 ;
    end
    Ratio_CC(ti) = Ratio_CC(ti) + CC/ ti ;



 gh = gh + 1 ;
    Conect2ty(gh, :, :) = Conect2ty(gh-1, :, :) ;


    for j = 1 : Size
        if j < m
            Conect2ty(gh, m, j) = P_Connection(m , j) ;
        end
        if j > m
            Conect2ty(gh, m, j) = P_Connection(m , j-1) ;
        end
    end

    for j = 1 : Size
        if j < n
            Conect2ty(gh, n, j) = P_Connection(n , j) ;
        end
        if j > n
            Conect2ty(gh, n, j) = P_Connection(n , j-1) ;
        end
    end


       
                Wealth(ti) =   Wealth(ti )  +  ( Out2( m ) + Out2( n ) )/(2 * Delay(ti) ) ; 
       
                AccumAllWealth(ti , m) = AccumAllWealth(ti-1 , m)  +  Out2( m ) ;
  AccumAllWealth(ti , n) = AccumAllWealth(ti-1 , n)  +  Out2( n ) ;

       for hh77 = 1 : Size
           if hh77 ~= m  && hh77 ~= n
           AccumAllWealth(ti , hh77) = AccumAllWealth(ti-1 , hh77)  +  0 ;
           end
       end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Infection
if  ti > TimeVirus 
  
                  %%% first agents update infaction then meet eachother:
                                                                                        %%%% ONLY CONNECTION HAPPENED AT ti  m and n
    if  Vaccinated(m) == 0 && Infected(n) == 1  %%%% Died(m) ==0  && 
        
        r = rand ;
        if r < InfectionChance
        Infected(m) = 1 ;                                                       %%% CAN BE PROBEBLISTIC
        end
    end
    
    if  Vaccinated(n) == 0 && Infected(m) == 1 %%%% Died(n) ==0  &&
        
        r = rand ;
        if r < InfectionChance
        Infected(n) = 1 ;
        end
        
    end
 
                                                        %%%%% Virus time
           %%%%% ISOLAED or DEAD agent COULD NoT PLAY
                                 %%% at each itteraton all the infected and isolated periods increase 1                                                 
    for yy11 = 1: Size
     if Infected(yy11) ==1   
        TimeInfected(yy11) = TimeInfected(yy11)  + 1 ;
    end
    end
    
    %%%%%%%%% AGENTS DECIDED TO BE ISOLATED ALTHOUGH WERE NOT ASKED TO PLAY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
     for yy11 = 1: Size
     if Isolated(yy11) ==1   
        TimeIsolated(yy11) = TimeIsolated(yy11)  + 1 ;
    end
    end
                                                                                   
                          %%%%   but might die or get isolated after two
                          %%%%   previous updates
        
          
%%%%%%%%%%%%%%  OUT of ISoLatioN
  
  
     for yy11 = 1: Size
 if TimeIsolated(yy11) >=  1  *     IsolationPeriod  
    TimeIsolated(yy11) = 0 ;  %%% in some cases might get new virus!!!
    Isolated(yy11) = 0  ;
end 
     end
    
                                     %%%%%%%%%%% TOP DOWN OUT of ISolAtiOn
     if  ti == TimeOuTiSolaTion
for oo11 = 1 : NumberIsolated
       Isolated(out0(oo11)) = 0 ;
end
     end
    
       
                                          %%%%%%%%%%%% Some agents DIE!
    for yy11 = 1: Size
                                                if    TimeInfected(yy11) > LifeTimeVirus   %%%%% rand <   (  exp(TimeInfected(yy11) - LifeTimeVirus )  )   %%%%   
                                                    
                                                    TimeInfected(yy11) = 0 ; 
                                                    TimeIsolated(yy11)= 0 ; 
                                                    
                                                     r = rand ;

                                                     if r < DeadChance
                                                         Died(yy11) = 1 ;
                                                         
                                                         Vaccinated(yy11) = 0 ;
                                                         Infected(yy11) = 0 ;
                                                         Isolated(yy11) = 0 ;

                                                     else
                                                         Died(yy11) = 0 ;
                                                         
                                                         Vaccinated(yy11) = 1 ;
                                                         Infected(yy11) = 0 ;
                                                         Isolated(yy11)= 0 ;
                                                     end
                                                end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END Infection
% 
% % 
% %                       %%%%%%%%%%%%%%%%%%   update Isolation FROM DESEASE    
%     
              %%%%%%% If n detects Symptoms in partner m, not dead then might go  to isolation

               if  Infected(m) == 1 &&  ( Isolated(n) == 0 && Died(n) == 0 )  &&  TimeInfected(m) >  SymptomPeriod 
                   
                   r = rand ;
                   if r < OthersSympt       
                    Isolated( n ) = 1 ;
                   end
               end
                      
               
               if  Infected(n) == 1 &&  ( Isolated(m) == 0 && Died(m) == 0 )   &&  TimeInfected(n) >  SymptomPeriod 
                   
                    r = rand ;
                   if r < OthersSympt        
                    Isolated( m ) = 1 ;
                   end
               end    
                      
               
               %%%%%%% Detencting symptoms in itself and might go to isolation   
               
              if  Infected(m) == 1 &&    ( Isolated(m) == 0 && Died(m) == 0 )  &&  TimeInfected(m) >  SymptomPeriod   %% &&  rand < TimeInfected(m) / LifeTimeVirus  %%%% rand <   (  exp(TimeInfected(m) - LifeTimeVirus )  )            
                   
                   r = rand ;
                   if r < SelfSympt      
                    Isolated( m ) = 1 ;
                   end
               end
                      
               
               if Infected(n) == 1 &&    ( Isolated(n) == 0 && Died(n) == 0 ) &&  TimeInfected(n) >  SymptomPeriod  %%&& rand < TimeInfected(n) / LifeTimeVirus  %%%  rand <   (  exp(TimeInfected(n) - LifeTimeVirus )  )            
                   
                    r = rand ;
                   if r < SelfSympt        
                    Isolated( n ) = 1 ;
                   end
               end    
                     
    Meann = 0 ;
 
    for lll = 1 : Size
        
        if Died(lll) == 0  && Isolated(lll) == 0
            if D_C_decision(lll)  == 2
        Meann = Meann + 1  ;
            end
        end
    end
       
    Meanfield(ti ) = Meanfield(ti ) + Meann / Size ;

if  D_C_decision(m) == 1  &&  D_C_decision(n) == 1
    C1C1 = C1C1 + 1 ;
end

MeanCC(ti) = MeanCC(ti) + C1C1 ;
    MeannDead = 0 ;
    for lll = 1 : Size
        MeannDead = MeannDead + Died(lll)  ;
    end
       
    MeanDied(ti ) = MeanDied(ti ) + MeannDead / Size ;

      MeannVac = 0 ;
    for lll = 1 : Size
        MeannVac = MeannVac + Vaccinated(lll)  ;
    end
       
    MeanVaccinated(ti ) = MeanVaccinated(ti ) + MeannVac / Size ;  
 
        MeannIsolInfe = 0 ;
        MeannIsolHealthy = 0 ;
    for lll = 1 : Size
        
             if Died(lll) == 0  
        if Isolated(lll) == 1
            
        if Infected(lll) == 1
        MeannIsolInfe = MeannIsolInfe + 1  ;
        else
            
     MeannIsolHealthy = MeannIsolHealthy + 1  ;

        end
        end
             end
    end
       
    MeannIsolatedInfec(ti ) = MeannIsolatedInfec(ti ) + MeannIsolInfe / Size ;
    MeannIsolatedHealthy(ti ) = MeannIsolatedHealthy(ti ) + MeannIsolHealthy / Size ;


              MeannFREEInfe = 0 ;
              MeannFREEHealth = 0 ;
    for lll = 1 : Size
       
      if Died(lll) == 0  
        if Isolated(lll) == 0
            
        if Infected(lll) == 1
        MeannFREEInfe = MeannFREEInfe + 1  ;
        
        else
     MeannFREEHealth = MeannFREEHealth + 1  ;
        end
        end
       
      end 
        
    end
       
    MeannFREEInfect(ti ) = MeannFREEInfect(ti ) + MeannFREEInfe / Size ;
        MeannFREEHealthy(ti ) = MeannFREEHealthy(ti ) + MeannFREEHealth / Size ;
  
        MeannInf = 0 ;
    for lll = 1 : Size
        MeannInf = MeannInf + Infected(lll)  ;
    end
       
    MeannInfected(ti ) = MeannInfected(ti ) + MeannInf / Size ;

end   %%%%% end time
 
end   % end ensemble


%  Demo
for ty = DemoTime : DemoTime

    weights(:, :) = Conect2ty(ty, :, :) ; % The intensity of the lines represents the propensity of the agents to connect/play with one another

    for tn  = 1 : Size
        Summ = cumsum(weights(tn, :)) ;
        weights(tn, :) = weights(tn, :) / Summ(Size) ;
    end

    G = digraph(weights) ;
    LWidths = 1 *G.Edges.Weight ;

    plot(G,'EdgeLabel',G.Edges.Weight,'LineWidth',LWidths)


    intensityValue = LWidths  ;
    OOOO =   intensityValue  ;
    hhh = plot(G,'LineWidth',LWidths) ;
    set(gcf,'color','w') ;
    cccc = 0 ;
    for uuu = 1 : Size
        uuu
        axis off
        for vvv = 1 : Size
            if uuu ~= vvv
                cccc =  cccc + 1 ;


                if   Infected(vvv) == 1 && Isolated(vvv) == 1
                    Col = [0.6 0.6 0] ; % infected and isolated, yellow dark
                
                elseif   Infected(vvv) == 1 && Isolated(vvv) == 0
                    Col = [0.9 0.9 0] ; % infected and not isolated, yellow light
                
                elseif Vaccinated(vvv) == 1
                    Col = [0 1 0] ; % Vaccinated green
                
                elseif Died(vvv) == 1
                    Col = [0.3 0.3 0.3] ; % Died black
                    
                else
                    if  Isolated(vvv) == 1
                    Col = [0 0.7 0.7] ; % Susceptible cyan light
                    elseif  Isolated(vvv) == 0
                    Col = [0 1 1] ; % Susceptible cyan
                    end
                end

                highlight(hhh, uuu, vvv, 'EdgeColor', [ (1 -OOOO(cccc, 1))  (1-OOOO(cccc, 1))  (1-OOOO(cccc, 1)) ], 'LineWidth', 2.5, 'MarkerSize', 8)  ;
               if D_C_decision(vvv) == 2
                highlight(hhh,  vvv, 'NodeColor', Col, 'Marker','o')  ;
               else
                highlight(hhh,  vvv,'NodeColor', Col, 'Marker','^','MarkerSize',8)  ;

               end
            end
        end
    end



    frame = getframe(gcf) ;
    writeVideo(mov,frame)
 
end
hold off

close(mov)


    MeannIsolatedInfec = MeannIsolatedInfec/ ENSRatio  ;
    MeannIsolatedHealthy =MeannIsolatedHealthy/ ENSRatio  ;
    MeannFREEInfect =MeannFREEInfect/ ENSRatio  ;
    MeannFREEHealthy =MeannFREEHealthy/ ENSRatio  ;
  
       
MeannIsolated = MeannIsolatedInfec + MeannIsolatedHealthy ;
                                                               FREEAll =  MeannFREEInfect +  MeannFREEHealthy ;             

Meanfield = Meanfield / ENSRatio ;
MeanCC = MeanCC  / ENSRatio  ;
MeanCDDC = MeanCDDC  / ENSRatio  ;
MeanDD = MeanDD  / ENSRatio  ;

MeanVaccinated = MeanVaccinated/ ENSRatio  ;
MeannIsolated = MeannIsolated/ ENSRatio  ;
MeanDied = MeanDied/ ENSRatio  ;
MeannInfected = MeannInfected/ ENSRatio  ;
Wealth = Wealth / ENSRatio  ;

AAAA =  (1.2/20) * Trial ;
    
TREGHV = zeros(Trial , 1) ;
TREGHV(:,1) = (1 + Meanfield(:,1))/2  ;

    ShortMeanfield = zeros(AAAA, 1) ;
    for ii = 1 : AAAA
        ShortMeanfield(ii) =                       0.1  *    TREGHV( ii*1  + TimeVirus  -  0.1*AAAA ) ;
    end
    
        ShortInfected = zeros(AAAA, 1) ;
    for ii = 1 : AAAA
        ShortInfected(ii) = MeannInfected( ii*1  + TimeVirus  -  0.1*AAAA ) ;
    end
    
            ShortIsolated = zeros(AAAA, 1) ;
    for ii = 1 : AAAA
        ShortIsolated(ii) = MeannIsolated( ii*1  + TimeVirus  -  0.1*AAAA ) ;
    end
    
               ShortVaccinated = zeros(AAAA, 1) ;
    for ii = 1 : AAAA
        ShortVaccinated(ii) = MeanVaccinated( ii*1  + TimeVirus  -  0.1*AAAA ) ;
    end
    
                  ShortDied = zeros(AAAA, 1) ;
    for ii = 1 : AAAA
        ShortDied(ii) = MeanDied( ii*1  + TimeVirus  -  0.1*AAAA ) ;
    end
    
    % plot(MeanDied,'DisplayName','MeanDied');hold on;plot(MeannInfected,'DisplayName','MeannInfected');plot(MeannIsolated,'DisplayName','MeannIsolated');plot(MeanVaccinated,'DisplayName','MeanVaccinated');hold off;
% plot(ShortDied,'DisplayName','ShortDied');hold on;plot(ShortInfected,'DisplayName','ShortInfected');plot(ShortIsolated,'DisplayName','ShortIsolated');plot(ShortVaccinated,'DisplayName','ShortVaccinated');hold off;

toc



% Decision making
function D = Decision(Prop)

g = length(Prop) ;
CS =   cumsum(Prop) ;

r = rand ;
for i = 1 : g
    if r < CS(i)
        D = i ;
        break
    end
end
end

% Payoff of the Prisoner's Dilemma game
function Payoff = Payoff_PD(DC1, DC2, s, p, tc)

if DC2 == 2
    NC = 1 ;
else
    NC = 0 ;
end

if DC1 == 2
    Payoff = NC * (1) + (1 - NC) * (-s) ;
else
    Payoff = NC * (1 + tc) + (1 - NC) * (p) ;
end

end


% Updating the propensities

function P = SA_Update(MechanismUsed, pp, Decision, Pay, Paybefore, Chi)

% This function updates the propensities of the decision mechanisms which
% contributed to the final decision of the agent
% MechanismUsed: 2 if the decision mechanism is used and 1 otherwise
% pp: The input propensities
% Decision: Indicates the index of the propensity that resulted in the decision
% Pay: Payoff of the agent at the current trial
% Paybefore: Payoff of the agent at the previous trial
% Chi: Constants relating payoff to the change of the propensity/probability

if MechanismUsed == 2 % If the decision mechanism has used (2) in the final decision making or not (1)
    g = length(pp) ;
    P = zeros(1 , g) ;

    Delta =  Chi * (Pay - Paybefore) / (abs(Pay) + abs(Paybefore)) ;
    if   Pay == Paybefore
        P = pp ;
    else

        P(Decision) = pp(Decision) + Delta ;

        for i = 1 : g
            if i ~= Decision
                P(i) = pp(i) - Delta/(g-1) ;
            end
        end

        % Boundary condition
        for j = 1 : g
            if  P(j) <= 0
                P(j) =  1e-10 ;
            end
            if  P(j) >= 1
                P(j) =  1-1e-10 ;
            end
        end

        % Normalization
        R = cumsum(P);
        P = P ./ R(g) ;

    end

else
    P = pp ;
end

end





