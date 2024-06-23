tic
% Korosh Mahmoodi 062224
clc;
clear all ;
close all ;

Trial = 2e6 ; % 2e6 ; 
                                            XR = 10 ; % for indexing number of initial isolated agents

MAXXX = 1000 ;   %%%%%% max value of the tendencies

delTrust = 0.3 * MAXXX ; % ratio of change in the trust tendency
delCD = 0.3 * MAXXX ;
delConnect = 0.3 * MAXXX ;

% parameters of the Prisoner's Dilemma game
s = 0 ;   
p = 0 ;
tt = 0.9 ; % temptation to cheat, tt < s + 1

                                        ENSRatio = 50 ; % number of simulations

TrialPayoff = zeros(XR, ENSRatio ) ; % payoff at last trial
TrialMeanf = zeros(XR, ENSRatio ) ; % mean field at last trial
TrialDead = zeros(XR, ENSRatio ) ; % number of died agents at last trial

Size  = 50 ; % number of the agents

InitialInfected = 20 ;  % number of initial infected agents

TimeVirus = 0.5 * Trial  ; % time when infection starts
LifeTimeVirus = 1e5 ; % each infected agent either dies or gets immuned after this period of time

TimeShotDown = TimeVirus ; % start of the top-down isolation

IsolationPeriod =  1.5 * LifeTimeVirus ; % time period that the agent goes to isolatoin

TimeOuTiSolaTion = TimeShotDown + IsolationPeriod ;  % Trial +1 ; %%% TimeShotDown + 0.5 * LifeTimeVirus ; %%% TimeShotDown   +  2 * LifeTimeVirus  ;  %%%

DeadChance = 0.4 ; % chance that infected agent dies after being infected for LifeTimeVirus
InfectionChance = 0.1 ; % chance that infection spreads to the pared agent

MeanDied = zeros( Trial, XR) ; % mean value of died agents
MeannInfected = zeros(Trial, XR) ;
MeanVaccinated = zeros(Trial, XR) ;
 
MeannIsolatedInfec = zeros(Trial, XR) ;
MeannIsolatedHealthy = zeros(Trial, XR) ;
MeannFREEInfect = zeros(Trial, XR) ;
MeannFREEHealthy = zeros(Trial, XR) ;
Meanfield = zeros(Trial, XR) ;

AccumPay = zeros(Trial, XR ) ; % accumulative sum of meanpayoff/delaytime
NumIso = zeros(XR, 1) ; % number of initial isolated agents

        Delay = zeros( Trial, XR) ; % number of tries in pairing process

for bv = 1 :     XR    
    CC = 0 ;
    bv
   NumIso(bv) = floor( ((bv-1)*0.1*Size) ) ;
    NumberIsolated = NumIso(bv) ;

    for zz = 1 : ENSRatio
        % zz
        AccumPayoff =  0 * ones(Size, 1) ; 

        Died = zeros( Size , 1) ;
        Infected = zeros( Size , 1) ;
        Vaccinated = zeros( Size , 1) ;
        Isolated = zeros( Size , 1) ;

        TimeInfected = zeros( Size , 1) ;
        TimeIsolated = zeros( Size , 1) ;

        %%%%%%%%%%%% for Selfish Algorithm

        CD_decision = zeros( Size , 1) ; %  Cooperation (C, or 1) or Defection decision (D, or -1)

        Pi1 = zeros( Size , Size) ; % tendency for decision C
        Qi1 = zeros( Size , Size) ; % tendency for decision D

        Pi2 = zeros( Size , Size) ; % updated tendency
        Qi2 = zeros( Size , Size) ;

        Out1 = zeros( Size , 1) ; % last payoff 
        Out2 = zeros( Size , 1) ; % current payoff

        Trust1 = zeros( Size , Size) ; % tendency for trust
        NotTrust1 = zeros( Size , Size) ; % tendency for not-trust

        Trust2 = zeros( Size , Size) ; % updated tendency
        NotTrust2 = zeros( Size , Size) ;

        Chance1 =  zeros(Size, Size); % tendency for pairing
        Chance2 =  zeros(Size, Size); % updated tendency

                                        % initial conditions
        for jjj = 1 : Size

            r = rand ;
            if r < 0.5
                CD_decision(jjj) = -1 ;
            else
                CD_decision(jjj) = -1 ;
            end
            for iii = 1 : Size
                if iii ~= jjj

                    Pi1(jjj , iii) = 0.1 * MAXXX ;
                    Qi1(jjj, iii)  = 0.9 * MAXXX ;

                    Trust1(jjj, iii) = 0  ;
                    NotTrust1(jjj, iii) =  MAXXX ;

                end
            end
        end


        for kkkk = 1 : Size

            r = rand ;
            if r < 0.5
                if CD_decision(kkkk) == 1
                    Out1( kkkk ) = 1 ;
                else
                    Out1( kkkk ) = -s ;
                end
            else

                if CD_decision(kkkk) == -1
                    Out1( kkkk ) = p ;
                else
                    Out1( kkkk ) = 1 + tt ;
                end
            end

        end


        for ia = 1 : Size
            for ib = 1 : Size
                if ia ~= ib
                    Chance1(ia, ib) =  MAXXX /2   ;
                end
            end
        end


        Meann = 0 ;
        for lll = 1 : Size
            if Died(lll) == 0 && Isolated(lll) == 0
            Meann = Meann + CD_decision(lll)  ;
            end
        end
        Meanfield(1, bv) = Meann / Size ;


        C1C1 = 0 ;
        gh =  1 ;
        for ti = 2 : Trial

           % initial infection
            if ti == TimeVirus  % a ratio randomly gets infected

                L0 = InitialInfected ;
                out = randperm(Size,L0)  ;

                for oo91 = 1 : InitialInfected
                    Infected(out(oo91)) = 1 ;
                end
            end

            tiDelay = 1 ;
                          %  top-down islolation
            if ti == TimeShotDown  
                n0 = Size  ;
                L0 = NumberIsolated ;
                out0 = randperm(n0,L0) ;

                for oo11 = 1 : NumberIsolated
                    Isolated(out0(oo11)) = 1 ;
                end
            end

           % picking agents to play
            m = floor( 1 + Size * rand) ;
            n = floor( 1 + Size * rand) ;

            while m == n
                n = floor( 1 + Size * rand) ;
            end

            Suum1 = 0 ;
            for to = 1 : Size
                Suum1 = Suum1 + Chance1( m , to  ) ;
            end

            P1 = Chance1( m , n ) /  Suum1 ; % propensity/probability of agent m to play with agent n

            Suum2 = 0 ;
            for to = 1 : Size
                Suum2 = Suum2 + Chance1( n , to  ) ;
            end

            P2 = Chance1( n , m ) /  Suum2 ;

            r1 = rand ;
            r2 = rand ;

            while r1 > P1  || r2 > P2          || Isolated(m) == 1  || Died(m) == 1 || Isolated(n)==1 || Died(n)==1

                                % update time duration of infection for agens
                if Infected(m) == 1
                    TimeInfected(m) = TimeInfected(m)  + 1 ;
                end

                % if Isolated(m) == 1
                %     TimeIsolated(m) = TimeIsolated(m)  + 1 ;
                % end

                if Infected(n) == 1
                    TimeInfected(n) = TimeInfected(n)  + 1 ;
                end

                % if Isolated(n) == 1
                %     TimeIsolated(n) = TimeIsolated(n)  + 1 ;
                % end

                m = floor( 1 + Size * rand) ;
                n = floor( 1 + Size * rand) ;

                while m == n
                    n = floor( 1 + Size * rand) ;
                end

                tiDelay =  tiDelay + 1 ;

                Suum1 = 0 ;
                for to = 1 : Size
                    Suum1 = Suum1 + Chance1( m , to ) ;   
                end
                P1 = Chance1( m , n ) / Suum1 ;

                Suum2 = 0 ;
                for to = 1 : Size
                    Suum2 = Suum2 + Chance1( n , to ) ;  
                end
                P2 = Chance1( n , m ) / Suum2 ;

                r1 = rand ;
                r2 = rand ;

            end

            Delay(ti) = tiDelay  ;  % number of attempts to pair

            if Infected(m) == 1
                TimeInfected(m) = TimeInfected(m)  + 1 ;
            end

            % if Isolated(m) == 1
            %     TimeIsolated(m) = TimeIsolated(m)  + 1 ;
            % end

            if Infected(n) == 1
                TimeInfected(n) = TimeInfected(n)  + 1 ;
            end

            % if Isolated(n) == 1
            %     TimeIsolated(n) = TimeIsolated(n)  + 1 ;
            % end


            % Cooperation or Defection decision

            Pmm = Pi1( m , n )/(Pi1(m , n) + Qi1(m , n)) ; % propensity for agent m to play as C with agend n
            Pnn = Pi1( n, m )/(Pi1(n , m) + Qi1(n, m )) ;


            r = rand ;
            if r < Pmm
                CD_decision(m) = 1 ;
            else
                CD_decision(m) = -1 ;
            end

            r = rand ;
            if r < Pnn
                CD_decision(n) = 1 ;
            else
                CD_decision(n) = -1 ;
            end

           
            %  imitation/Trust decision

            Votemm = CD_decision(m) ;
            Votenn = CD_decision(n) ;

            imitnn = 0 ;

            R0000 =  Trust1( n , m) / ( Trust1( n, m) +  NotTrust1( n, m) ) ; % propensity for agent n to trust the decision of agent m

            r = rand ;
            if r < R0000
                imitnn = 1 ;

                CD_decision(n) =   Votemm ;

            end

            Suum2 = 0 ;
            for to = 1 : Size
                Suum2 = Suum2 + Chance1( 1 , to  ) ;   
            end

            imitmm = 0 ;

            R0000 =  Trust1( m , n)  /  ( Trust1( m, n) +  NotTrust1( m, n) )  ;



            r = rand ;
            if r < R0000
                imitmm = 1 ;

                CD_decision(m ) =   Votenn ;


            end


           % evaluating the payoff

            if CD_decision(m) == 1
                ggnn = 1 ;
            else
                ggnn = 0 ;
            end

            if CD_decision(n) == 1
                ggmm = 1 ;
            else
                ggmm = 0 ;
            end

            if CD_decision(m) == 1
                Out2( m ) = ggmm * (1) + (1 - ggmm) * (-s);
            else
                Out2( m ) = ggmm * (1 + tt) + (1 - ggmm) * (p);
            end


            if CD_decision(n) == 1
                Out2( n ) = ggnn * (1) + (1 - ggnn) * (-s);
            else
                Out2( n ) = ggnn * (1 + tt) + (1 - ggnn) * (p);
            end

           %   updating the tendencies of Trust SAT

            if imitnn == 1

                if    Out2( n) == Out1( n)
                    Trust2( n , m) = Trust1(n, m) ;
                    NotTrust2( n, m) = NotTrust1(n, m) ;
                else
                    if  Out2( n) > Out1( n)
                        Trust2( n, m) =  Trust1( n, m) + delTrust * abs(Out2( n) - Out1( n)) / (abs(Out2( n)) + abs(Out1( n)) );
                        NotTrust2( n, m) = NotTrust1( n, m ) - delTrust * abs(Out2( n) - Out1( n)) / (abs(Out2( n)) + abs(Out1( n)) );
                    else
                        Trust2( n, m) =  Trust1( n, m) - delTrust * abs(Out2( n) - Out1( n)) / (abs(Out2( n)) + abs(Out1( n)) );
                        NotTrust2( n, m) = NotTrust1( n, m ) + delTrust * abs(Out2( n) - Out1( n)) / (abs(Out2( n)) + abs(Out1( n)) );
                    end

                end
            else

                if    Out2( n) == Out1( n)
                    Trust2( n , m) = Trust1(n, m) ;
                    NotTrust2( n, m) = NotTrust1(n, m) ;
                else

                    if  Out2( n) > Out1( n)
                        Trust2( n, m) = Trust1( n, m ) - delTrust * abs(Out2( n) - Out1( n)) / (abs(Out2( n)) + abs(Out1( n)) );
                        NotTrust2( n, m) =  NotTrust1( n, m ) + delTrust    * abs(Out2( n) - Out1( n)) / (abs(Out2( n)) + abs(Out1( n)) );
                    else
                        Trust2( n, m) = Trust1( n, m ) + delTrust * abs(Out2( n) - Out1( n)) / (abs(Out2( n)) + abs(Out1( n)) );
                        NotTrust2( n, m) =  NotTrust1( n, m ) - delTrust    * abs(Out2( n) - Out1( n)) / (abs(Out2( n)) + abs(Out1( n)) );
                    end

                end
            end


            if Trust2( n, m) < 0   || NotTrust2( n, m) > MAXXX
                Trust2( n, m) = 0 ;
                NotTrust2( n, m) = MAXXX ;
            end

            if NotTrust2( n, m) < 0  ||  Trust2( n, m) > MAXXX
                NotTrust2( n, m) = 0 ;
                Trust2( n, m) = MAXXX ;
            end

            if imitmm == 1
                if    Out2( m) == Out1( m)
                    Trust2( m , n) = Trust1(m, n) ;
                    NotTrust2( m, n) = NotTrust1(m, n) ;
                else

                    if  Out2( m) > Out1( m)
                        NotTrust1( m, n) =  Trust1( m, n) + delTrust   * abs(Out2( m) - Out1( m)) / (abs(Out2( m)) + abs(Out1( m)) );
                        NotTrust2( m, n) = NotTrust1( m, n )  - delTrust * abs(Out2( m) - Out1( m)) / (abs(Out2( m)) + abs(Out1( m)) );
                    else
                        Trust2( m, n) =  Trust1( m, n) - delTrust     * abs(Out2( m) - Out1( m)) / (abs(Out2( m)) + abs(Out1( m)) );
                        NotTrust2( m, n) = NotTrust1( m, n )  + delTrust * abs(Out2( m) - Out1( m)) / (abs(Out2( m)) + abs(Out1( m)) );
                    end

                end
            else


                if    Out2( m) == Out1( m)
                    Trust2( m , n) = Trust1(m, n) ;
                    NotTrust2( m, n) = NotTrust1(m, n) ;
                else

                    if  Out2( m) > Out1( m)
                        Trust2( m, n) = Trust1( m, n ) - delTrust  * abs(Out2( m) - Out1( m)) / (abs(Out2( m)) + abs(Out1( m)) );
                        NotTrust2( m, n) =  NotTrust1( m, n ) + delTrust  * abs(Out2( m) - Out1( m)) / (abs(Out2( m)) + abs(Out1( m)) );
                    else

                        Trust2( m, n) = Trust1( m, n ) + delTrust  * abs(Out2( m) - Out1( m)) / (abs(Out2( m)) + abs(Out1( m)) );
                        NotTrust2( m, n) =  NotTrust1( m, n ) - delTrust * abs(Out2( m) - Out1( m)) / (abs(Out2( m)) + abs(Out1( m)) );
                    end
                end
            end


            if Trust2( m, n) < 0  ||   NotTrust2( m, n) > MAXXX
                Trust2( m, n) = 0 ;
                NotTrust2( m, n) = MAXXX ;
            end

            if NotTrust2( m, n) < 0  ||  Trust2( m, n) > MAXXX
                NotTrust2( m, n) = 0 ;
                Trust2( m, n) = MAXXX ;
            end

           %   updating the tendencies of C or D decision SAL

            if imitmm == 0

                if CD_decision(m) == 1

                    T1 =  Out2( m)  ;
                    OUT0 =  Out1( m)  ;

                    if   T1 ==   OUT0
                        Pi2( m , n) =  Pi1(m, n)  ;
                        Qi2( m, n) =  Qi1( m, n) ;
                    else

                        if   T1 >   OUT0
                            Pi2( m , n) =  Pi1(m, n)  +    delCD   * abs(Out2( m) - Out1( m)) / (abs(Out2( m)) + abs(Out1( m)) );  %%%%   + 1 ;
                            Qi2( m, n) =  Qi1( m, n)  - delCD   * abs(Out2( m) - Out1( m)) / (abs(Out2( m)) + abs(Out1( m)) );
                        else
                            Pi2( m , n) =  Pi1(m, n)  -    delCD   * abs(Out2( m) - Out1( m)) / (abs(Out2( m)) + abs(Out1( m)) );  %%%%   + 1 ;
                            Qi2( m, n) =  Qi1( m, n)  + delCD   * abs(Out2( m) - Out1( m)) / (abs(Out2( m)) + abs(Out1( m)) );
                        end
                    end
                end



                if CD_decision(m) == -1
                    T1 =  Out2( m)  ;
                    OUT0 =  Out1( m)  ;

                    if   T1 ==   OUT0
                        Qi2( m , n) = Qi1( m , n)  ;
                        Pi2(m, n) =  Pi1(m, n)     ;
                    else

                        if   T1 >   OUT0
                            Qi2( m , n) =  Qi1( m , n) + delCD   * abs(Out2( m) - Out1( m)) / (abs(Out2( m)) + abs(Out1( m)) );
                            Pi2(m, n) =  Pi1(m, n)  - delCD  * abs(Out2( m) - Out1( m)) / (abs(Out2( m)) + abs(Out1( m)) )   ;  %%%%   + 1 ;
                        else
                            Qi2( m , n) =  Qi1( m , n) - delCD   * abs(Out2( m) - Out1( m)) / (abs(Out2( m)) + abs(Out1( m)) )  ;
                            Pi2(m, n) =  Pi1(m, n)  + delCD    * abs(Out2( m) - Out1( m)) / (abs(Out2( m)) + abs(Out1( m)) )  ;  %%%%   + 1 ;
                        end
                    end
                end

            else
                Qi2( m , n) = Qi1( m , n)  ;
                Pi2(m, n) =  Pi1(m, n)     ;  
            end


            if  Pi2( m , n)  < 0  ||  Qi2( m , n)  > MAXXX
                Pi2( m , n)  = 0   ;
                Qi2( m , n)  = MAXXX   ;
            end

            if  Qi2( m , n)  < 0 ||  Pi2( m , n)  > MAXXX
                Qi2( m , n)  = 0   ;
                Pi2( m , n)  = MAXXX   ;
            end



            if imitnn == 0

                if CD_decision(n) == 1
                    T1 =  Out2(n)  ;
                    OUT0 =  Out1(n)  ;

                    if   T1 ==   OUT0
                        Pi2(n , m) =  Pi1(n, m)    ;
                        Qi2(n, m) = Qi1(n, m)  ;
                    else
                        if   T1 >   OUT0
                            Pi2(n , m) =  Pi1(n, m)  +    delCD   * abs(Out2( n) - Out1( n)) / (abs(Out2( n)) + abs(Out1( n)) );  %%%%   + 1 ;
                            Qi2(n, m) =  Qi1(n, m)  - delCD  * abs(Out2( n) - Out1( n)) / (abs(Out2( n)) + abs(Out1( n)) ) ;
                        else
                            Pi2(n , m) =  Pi1(n, m)  -    delCD  * abs(Out2( n) - Out1( n)) / (abs(Out2( n)) + abs(Out1( n)) ) ;  %%%%   + 1 ;
                            Qi2(n, m) =  Qi1(n, m)  + delCD   * abs(Out2( n) - Out1( n)) / (abs(Out2( n)) + abs(Out1( n)) ) ;
                        end
                    end
                end

                if CD_decision(n) == -1
                    T1 =  Out2(n)  ;
                    OUT0 =  Out1(n)  ;
                    if   T1 ==   OUT0
                        Qi2(n, m) = Qi1(n , m) ;
                        Pi2(n, m) =  Pi1(n, m)   ;
                    else
                        if   T1 >   OUT0
                            Qi2(n, m) =  Qi1(n , m) + delCD   * abs(Out2( n) - Out1( n)) / (abs(Out2( n)) + abs(Out1( n)) ) ;
                            Pi2(n, m) =  Pi1(n, m)   - delCD   * abs(Out2( n) - Out1( n)) / (abs(Out2( n)) + abs(Out1( n)) ) ;  %%%%   + 1 ;
                        else
                            Qi2(n, m) =  Qi1(n , m) - delCD   * abs(Out2( n) - Out1( n)) / (abs(Out2( n)) + abs(Out1( n)) ) ;
                            Pi2(n, m) =  Pi1(n, m)   + delCD  * abs(Out2( n) - Out1( n)) / (abs(Out2( n)) + abs(Out1( n)) ) ;  %%%%   + 1 ;
                        end

                    end

                end
            else
                Qi2(n, m) = Qi1(n , m) ;
                Pi2(n, m) =  Pi1(n, m)   ;  

            end
            if  Pi2(n, m)  < 0   ||  Qi2(n, m)  > MAXXX
                Pi2(n, m)  = 0   ;
                Qi2(n, m)  = MAXXX   ;
            end
            if  Qi2(n, m)  < 0 || Pi2(n, m)  > MAXXX
                Qi2(n, m)  = 0   ;
                Pi2(n, m)  = MAXXX   ;
            end

           %   updating the tendencies of to whom to play decision SAC

            if    Out2( m) == Out1( m )
                Chance2( m , n ) =  Chance1( m , n ) ;
            else

                if   Out2( m) > Out1( m )
                    Chance2( m , n ) = Chance1( m , n ) + delConnect    * abs(Out2( m) - Out1( m)) / (  abs(Out2( m)) + abs(Out1( m))) ;
                else
                    Chance2( m , n ) = Chance1( m , n ) - delConnect   * abs(Out2( m) - Out1( m)) / (  abs(Out2( m)) + abs(Out1( m)));
                end
            end

            if    Out2( n) == Out1( n)
                Chance2( n , m ) =  Chance1( n , m ) ;
            else

                if   Out2( n) > Out1( n)
                    Chance2( n , m ) = Chance1( n , m ) + delConnect    * abs(Out2( n) - Out1( n)) / (  abs(Out2( n)) + abs(Out1( n)));
                else
                    Chance2( n , m ) = Chance1( n , m ) - delConnect   * abs(Out2( n) - Out1( n)) / (  abs(Out2( n)) + abs(Out1( n)));
                end
            end

            if   Chance2( m , n ) <= 0
                Chance2( m ,n ) =  Chance1( m , n )  ;
            end

            if   Chance2( n , m ) <= 0
                Chance2( n , m ) =  Chance1( n , m ) ;
            end




                                % mean field
            Meann = 0 ;
            for lll = 1 : Size
                            if Died(lll) == 0 && Isolated(lll) == 0

                Meann = Meann + CD_decision(lll)  ;
                            end
            end
            Meanfield(ti, bv ) = Meann / Size ;


            % updating the current variables as past (for next trial)

            Chance1( m , n ) =  Chance2( m , n )  ;
            Chance1( n , m ) =  Chance2( n , m )   ;

            Pi1(m , n) =  Pi2(m , n) ;
            Qi1(m , n) =  Qi2(m , n) ;

            Out1(m) = Out2(m) ;

            Trust1( m , n) = Trust2( m , n) ;
            NotTrust1( m , n) = NotTrust2( m , n) ;

            Pi1(n , m) =  Pi2(n, m) ;
            Qi1(n, m) =  Qi2(n, m) ;

            Out1(n) = Out2(n) ;
            Trust1( n, m) = Trust2( n, m) ;
            NotTrust1( n, m) = NotTrust2( n, m) ;

            % accumulative payoff
        if ti >= TimeVirus  && ti <= TimeVirus +   1e5

            AccumPay(ti, bv) =   AccumPay(ti-1, bv)  +  ( Out2( m ) + Out2( n ) )/(2 * Delay(ti) ) ;
        end


                                              %  infection
            if  ti > TimeVirus

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

              

                                  %  out of isolation
                if  ti == TimeOuTiSolaTion
                    for oo11 = 1 : NumberIsolated
                        Isolated(out0(oo11)) = 0 ;
                    end
                end


                                             % agents die or immune
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



            MeannDead = 0 ;
            for lll = 1 : Size
                MeannDead = MeannDead + Died(lll)  ;
            end

            MeanDied(ti, bv ) =  MeannDead / Size ;

            MeannVac = 0 ;
            for lll = 1 : Size
                MeannVac = MeannVac + Vaccinated(lll)  ;
            end

            MeanVaccinated(ti, bv ) =  MeannVac / Size ;

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

            MeannIsolatedInfec(ti, bv ) =  MeannIsolInfe / Size ;
            MeannIsolatedHealthy(ti, bv ) = MeannIsolHealthy / Size ;


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

            MeannFREEInfect(ti, bv ) =  MeannFREEInfe / Size ;
            MeannFREEHealthy(ti, bv ) = MeannFREEHealth / Size ;

            MeannInf = 0 ;
            for lll = 1 : Size
                MeannInf = MeannInf + Infected(lll)  ;
            end

            MeannInfected(ti, bv ) =  MeannInf / Size ;


        end   %%%%% end time

        TrialDead(bv, zz) = MeanDied(Trial, bv) ;  % Meanfield(Trial, bv) ; %  TrialPayoff(bv) + AccumPay(Trial, bv) ;
        TrialMeanf(bv, zz) =  Meanfield(Trial, bv) ; %  TrialPayoff(bv) + AccumPay(Trial, bv) ;Ratio_CC

        TrialPayoff(bv, zz) = AccumPay(TimeVirus+1e5, bv) ;  %  AccumPay(Trial, bv) ;

    end   % end ensemble

end  % end ratio of number of isolated

MeannIsolatedInfec = MeannIsolatedInfec/ ENSRatio  ;
MeannIsolatedHealthy =MeannIsolatedHealthy/ ENSRatio  ;
MeannFREEInfect =MeannFREEInfect/ ENSRatio  ;
MeannFREEHealthy =MeannFREEHealthy/ ENSRatio  ;

MeannIsolated = MeannIsolatedInfec + MeannIsolatedHealthy ;

Meanfield = Meanfield / ENSRatio ;

MeanVaccinated = MeanVaccinated/ ENSRatio  ;
MeannIsolated = MeannIsolated/ ENSRatio  ;
MeanDied = MeanDied/ ENSRatio  ;
MeannInfected = MeannInfected/ ENSRatio  ;



Mfield=mean(Meanfield') ;
Mfield=Mfield' ;


st = TimeVirus ;
del = 2e4 ;
Mdied = mean(MeanDied') ;
Mdied = Mdied' ;
MshortDied = Mdied(st-1e3:st+del, 1) ;

Minfect = mean(MeannInfected') ;
Minfect = Minfect' ;
MshortInfec = Minfect(st-1e3:st+del, 1) ;

Miso = mean(MeannIsolated') ;
Miso = Miso' ;
MshortIso = Miso(st-1e3:st+del, 1) ;

Mvac = mean(MeanVaccinated') ;
Mvac = Mvac' ;
MshortVac = Mvac(st-1e3:st+del, 1) ;

range = st-1e3:st+del ;

Matrixx = zeros(length(MshortVac), 5) ;

Matrixx(:, 1) =  range; 
Matrixx(:, 2) =  MshortDied; 
Matrixx(:, 3) = MshortInfec ; 
Matrixx(:, 4) = MshortIso ; 
Matrixx(:, 5) = MshortVac ; 

plot(Mfield) ;

xx = mean(TrialDead') ; xx = xx' ;
yy = mean(TrialPayoff'); yy = yy' ;


toc