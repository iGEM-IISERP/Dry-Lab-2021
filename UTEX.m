% ------------------ UTEX 2973 SUCROSE OPTIMIZATION MODEL ----------------

% Initiating the package and assigning a good solver
% Run this block only the first time you run the code, comment is out for
% all the subsequent trials
%-----------------------------------------------------------------------
%{
global TUTORIAL_INIT_CB;
if ~isempty(TUTORIAL_INIT_CB) && TUTORIAL_INIT_CB==1
initCobraToolbox(false) % false, as we don't want to update
changeCobraSolver('gurobi','all');
end
%}

% -----------------------------------------------------------------------

%{
% Read the model
% Note: this is specific to my file. For your file to be read
% just change the path and file name accordingly.
model = readCbModel('UTEX.xml');


% Apparently the model did not have any reaction to produce sucrose 
% Therefore we had to add the SPS and SPP reactions from the PCC 7942 model
model = addReaction(model, 'SPS', 'reactionFormula',...
    'cpd00072[c] + cpd00026[c] -> cpd00014[c] + cpd00067[c] + cpd00077[c]'...
    , 'lowerBound', 0, 'upperBound', 1000);


% cpd00077 is sucrose 6 phosphate which is the compund from which sucrose
% is produced
model = addReaction(model, 'SPP', 'reactionFormula',...
    'cpd00001[c] + cpd00077[c] -> cpd00009[c] + cpd00076[c]', ...
    'lowerBound', 0, 'upperBound', 1000);


% sink reaction to account for sucrose accumulation within the cell to
% counter act osmotic pressure
model = addReaction(model, 'Sucrose_acc', 'metaboliteList', ...
    {'cpd00076[c]'}, 'stoichCoeffList', [-1], 'reversible', false);


% Sucrose is made in the cytoplasm ie [c] and therefore we need a
% transporter to transport it outside... heterologus expression of cscB is
% the answer
model = addReaction(model, 'cscB', 'reactionFormula',...
    'cpd00076[c] + cpd00067[c] -> cpd00076[p] + cpd00067[p]', ...
    'lowerBound', -1000, 'upperBound', 1000);


% Variables used for graphs
markers = ['d'; '^'; 'o'; 's'; 'p'];
markercolor = ['#1A85FF';'#28C79B';'#FCF462';'#E43995';'#ffffff'];


% -----------------------------------------------------------------------
% Please note that this code contains reaction numbers specific to
% the model used while creating the code, 890 - biomass and 1074 - sucrose
% export. Change these numbers to the appropriate ones in your model
% before usin gthis code 
% -----------------------------------------------------------------------
 


 

% IGNORE THIS PART IF YOU ARE NOT SIMULATING SALT STRESS. THIS CODE IS TO
% ACCOUNT FOR HOW THE BACTERIA BEHAVES WHEN PUT IN AN ENVIORNMENT THAT
% FORCES IT TO MAKE METABOLIC ADJUSTMENTS.
    
% ----------------------------- M O M A ----------------------------------

% Co2 uptake rates as mentioned in Alberthy et al 2017
model.lb(895) = -12.2;


% These variables store values of CO2 uptakes rates, cscB efficiency
% values, sucrose production vlaues, and growth rates in a way that makes
% generating graphs easier
co_2 = zeros(1, 400);
cscB = zeros(1, 400);
suc = zeros(1, 400);
bm = zeros(1, 400);
counter = 0;
indices = 1;



for k = 1:20
    % From literature the uptake value for CO2 exchange was found to be -12.2
    model.lb(895) = -14.2 + (2*(k/10));

    

    % Calculating the growth rates for the organism outside of salt stress
    sol = optimizeCbModel(model);


    
    % Store the original flux values for MoMA
    v_original_biomass = sol.v;



    % for a 150mM NaCl medium the internal sucrose conc is around 0.17mmols 
    % per gram of biomass 
    model = addRatioReaction(model, {'Biomass_Auto_2973', 'Sucrose_acc'},...
        [0.17 1]);

    

    for j = 1:20
        % cscB has shown a transport efficiency of 90% therefore the sucrose
        % transport is 9 times the sucrose accumulation. song et al
        % here we take values from 80% to 90%
        model = addRatioReaction(model, {'Sucrose_acc', 'EX_Sucrose'},...
            [(9.5 - (j/10)) 1]);



        % changing the objective function 
        model.c(890) = 0;
        model.c(1074) = 1;

    
        % an array to store the different flux differences due to change in
        % lowerbound of biomass production
        delta = 1000 * ones(1, 100);

    
        % Loop to set different values as lower bounds of biomass production
        % reaction and then perform FBA to optimize for sucrose and then to store
        % the flux difference from the base case, ie. optmized for growth without
        % taking in salt stress
        for i = 1:100
    
            % 0.3033 is the FBA solution when you maximize for biomass
            model.lb(890) = 0.3033 * (i/100);
    
            sol = optimizeCbModel(model);
    
            % To check if the solver was successful, if yes, then we store
            % all the values
            if sol.stat ~= 0
                a = (0.5 * sol.v' * sol.v);
                b = (v_original_biomass' * sol.v);
                delta(i) = a - b;
                counter = counter + 1;
            end
    
        end
        

        % Finding the min flux difference 
        [val , index] = min(delta);

        % Assigning the corresponding lb value to biomass reaction
        model.lb(890) = 0.3033 * (index/1000);

        fprintf('CO2 intake bound: %4.2f\n', -14.2 + (2*(k/10)))
        fprintf('cscB Efficiency: %2f percentage\n', 10 * (9.5 - (j/10))) 

        % Optimizing for product 
        sol = optimizeCbModel(model);

    
        % This is to get live results as the code is running, you can
        % comment this part out as it was mostly used for bug fixes
        fprintf('Sucrose excretion rate: %6.4f\n',sol.f)
        fprintf('Biomass Production rate: %6.4f\n\n\n\n', sol.v(890))
        
        
        % Storing values to generate graphs
        co_2(indices) = -14.2 + (2*(k/10));
        cscB(indices) = (9.5 - (j/10));
        suc (indices) = sol.f;
        bm (indices) = sol.v(890);
        indices = indices + 1;

        
        % To update constraints we remove the previous
        % New ones are added at the start of every loop
        model = removeCOBRAConstraints(model, ...
            {'Ratio_Sucrose_acc_EX_Sucrose'});

    end
    
    
     % To update constraints we remove the previous
     % New ones are added at the start of every loop
     model = removeCOBRAConstraints(model, ...
            {'Ratio_Biomass_Auto_2973_Sucrose_acc'});
        
        
end
%}

figure;


% To plot a graph that contains all the variables we have kept a tab on.
scatter3(co_2, cscB, bm, 50, suc, 'filled')
title('Sucrose Production [mmol/gDW/hr]', 'Color', 'white')
xlabel('CO_2 Uptake rate [mmol/gDW/hr]')
ylabel('cscB efficiency [per10]')
zlabel('Biomass Production [/hr]')
hcb = colorbar;
hcb.Color =  'white';
hcb.Title.String = "Sucrose Production";
hcb.Title.Color = 'white';


%{
counter = 1;
cscb = [94, 90, 85, 80, 75];
set(gcf,'position',[10,10,800,400])
for i=[1,4,9,14,19]
    for j=1:20
        y(j) = bm((j-1)*20 + i);
        x(j) = co_2((j-1)*20 + i);
    end
    plot(x, y, markers(counter), 'Marker', markers(counter),'Color', markercolor(counter,:), ...
        'MarkerFaceColor',markercolor(counter,:),...
        'MarkerSize',4, 'DisplayName',strcat("cscB efficiency is ",string(cscb(counter))))
   hold on
   counter = counter + 1;
end
%}

%xlabel('Co_2 uptake rate [mmol/gDW/hr]')
%ylabel('Biomass [hr^-^1]') 
%title('Biomass vs CO_2 uptake rates', 'Color', 'white')
ax = gca; % current axes
ax.FontSize = 12;
ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
ax.XAxis.Color='white';
ax.YAxis.Color='white';
ax.ZAxis.Color='white';
ax.GridColor = 'white';
ax.Color = '#212529';
ax.LineWidth = 2.0;
ax.AmbientLightColor = [1 0 1];
ax.Box = 'off';
set(gcf,'Color','#212529');
view(45,25)
hold off
%legend('location', 'southeastoutside', 'EdgeColor', '#fff', 'TextColor', '#fff');
exportgraphics(gcf, "3_Variables.jpg", 'BackgroundColor', '#212529');
savefig("3_Variables");


% ---------------------- E N D  O F  M O M A -----------------------------









%{
% --------------------------- F S E O F ----------------------------------

%                              READ ME
% To use the model without simulating salt stress comment out the two
% ratioReactions mentioned below. Furthur below you will see the code for
% OptKnock without salt stress. I personally have not tried running this
% code without salt stress so i'm not sure if this will give you any
% outputs. 



% Completing the model with constraints based on its salt stress
% Assume a cscB efficieny of 90%
%model = addRatioReaction(model, {'Sucrose_acc', 'EX_Sucrose'},...
    %[9 1]);


% for a 150mM NaCl medium the internal sucrose conc is around 0.17mmols 
% per gram of biomass 
%model = addRatioReaction(model, {'Biomass_Auto_2973', 'Sucrose_acc'},...
       % [0.17 1]);
    

% Highest uptake rate found in literature 
model.lb(895) = -12.2;


% To find the max possible sucrose production rate
model.c(890) = 0;
model.c(1074) = 1;
v_suc_star = optimizeCbModel(model).v(1074);


% To find max possible growth rate
model.c(890) = 1;
model.c(1074) = 0;
v_bio_star = optimizeCbModel(model).v(890);


% Variable to store all the flux distributions
flux = zeros(length(model.rxns()), 11);


% Variable to store reactions that can be over expressed
react = zeros(1, length(model.rxns()));
imp = zeros(1, length(model.rxns()));
counter = 1;

% We will be using this variable for all future purposes in this section of
% the code so we add the "WT" flux distribution into the variable so we
% won't have to work with two different variables
flux(:, 1) = optimizeCbModel(model).v;



% Setting various restrictions on the sucrose production rates and 
% maximising for biomass. 
for i = 1:10
   model.lb(1074) = 0.1 * i * v_suc_star;
   sol = optimizeCbModel(model);
   flux(:, i+1) = sol.v; 
end



% reactions that increase with increase in sucrose production rates
for i = 1:length(model.rxns())
   if  flux(i, 3:11) > flux(i,2)
       react(counter) = i;
       %imp(counter) = flux(i, 11) - flux(i, 2);
       imp(counter) = max(flux(i,2:11)) - flux(i,2);
       counter = counter + 1;
   end
end


% To arrange all the reactions in order of impact
% Sort Method - Bubble Sort
for i = 1:counter - 1
    for j = 1:counter - i - 1
        if imp(j) < imp(j+1)
            
            temp = imp(j);
            imp(j) = imp(j+1);
            imp(j+1) = temp;
            
            temp = react(j);
            react(j) = react(j+1);
            react(j+1) = temp;
        end
    end
end


% Plot graphs 
% 10 graphs are plotted. The top 10 modifications based on impact is
% plotted


nil = 0;

for i = 1:counter-1
    
    
    figure;
    
    
    gprs = string(findGPRFromRxns(model, model.rxns(react(i))));
    newstr = split(gprs);
    if gprs == ""
        nil = nil + 1;
    end
    
    model.rxnNames(react(i));
    gprs;
    fprintf('\n')
    
    for j = 1:length(newstr)
        if newstr(j) ~= "or" && newstr(j) ~= "and" && newstr(j) ~= "(" && ...
                newstr(j) ~= ")"
            gene(end + 1) = newstr(j);
        end
    end
    
    
    % The x axis is the direction of increase of sucrose production
    % Y axis is the value of the flux of that particular reaction for
    % different sucrose fluxes
    
    
    set(gcf,'position',[10,10,800,400])
    %plot([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1],...
        %flux(react(i),:),'o', 'color', '#28C79B', 'DisplayName',...
        %char(model.rxnNames(react(i))),'MarkerFaceColor','#28C79B',...
       % 'MarkerSize',8, 'Marker','o')
    
    
    scatter3(((1:11) * 0.1) - 0.1 ,flux(react(i),:),flux(890,:),50,...
        'filled');
    xlabel('Sucrose Values [fraction of maximum sucrose flux]')
    ylabel('Reaction Flux [mmol/gDW/hr]')
    zlabel('Biomass Values [hr^-^1]')
    title(char(model.rxnNames(react(i))),'Color', 'White')
    %title('Reaction Fluxes vs Fraction of Maximum Sucrose','Color','white') 
    %colorbar('Color', 'white');
    
    ax = gca; % current axes
    ax.FontSize = 12;
    ax.TickDir = 'out';
    ax.TickLength = [0.02 0.02];
    ax.XAxis.Color='white';
    ax.YAxis.Color='white';
    ax.ZAxis.Color='white';
    ax.GridColor = 'white';
    ax.Color = '#212529';
    ax.LineWidth = 2.0;
    ax.AmbientLightColor = [1 0 1];
    ax.Box = 'off';
    set(gcf,'Color','#212529');
   % exportgraphics(gcf, strcat(string(i),".jpg"),'BackgroundColor', '#212529')
    savefig(string(i))

    %savefig(string(i))
    
end

%legend('location', 'southeastoutside', 'EdgeColor', '#fff', 'TextColor', '#fff');
%exportgraphics(gcf, strcat(string(i),".jpg"),'BackgroundColor', '#212529')
%savefig(string(i))

length(gene(gene==""));
gene = gene(gene~="");
gene = gene(gene~="a");
ids = {'a'};
for k = 1:length(gene)
    fprintf('\n\n\n\n')
    gene(k);
    for j = 1:1180             
        gprs_test = findGPRFromRxns(model, model.rxns(j));
        if contains(gprs_test, gene(k))
            ids {end + 1} = char(model.rxnNames(j));
            model.rxns(j);
            model.rxnNames(j);
        end
    end
end

% ------------------------ E N D  O F  F S E O F -------------------------
%}






%{
% ------------------------ O P T K N O C K -------------------------------

%                              READ ME
% To use the model without simulating salt stress comment out the two
% ratioReactions mentioned below. Again, I havent ran this code without
% salt stress so i'm not sure if you'll get any non-zero values. 


% Completing the model with constraints based on its salt stress
% Assume a cscB efficieny of 90%
%model = addRatioReaction(model, {'Sucrose_acc', 'EX_Sucrose'},...
    %[9 1]);


% for a 150mM NaCl medium the internal sucrose conc is around 0.17mmols 
% per gram of biomass 
%model = addRatioReaction(model, {'Biomass_Auto_2973', 'Sucrose_acc'},...
        %[0.17 1]);


% Highest uptake rate found in literature 
model.lb(895) = -12.2;


    
% Finding the set of reactions to delete
% Finding all the reactions that contain sucrose in either product or
% reactant side Note:- cpd00076 is the MetID for sucrose
rxn_set = findRxnsFromMets(model, {'cpd00076[e]', 'cpd00076[p]',...
    'cpd00076[c]'});


% Finding all the metabolites from the list of reactions with sucrosoe in 
% them
%met_set = findMetsFromRxns(model, rxn_set); 


% Finding all the reactions with the metabolites from the previous list.
% These are basically all the reaction either directly or indirectly
% invovlved with sucrose
SelRxn = model.rxns;


% before optKnock we store the flux distribution of the "WT" 
model.c(890) = 1;
model.c(1074) = 0;
solWT = optimizeCbModel(model);




%{
% ------------------------ NO SALT STRESS --------------------------------



% This variable keeps record of the target reaction and the limit on the
% number of knockouts 

% Syntax - options = struct('targetRxn', {{sucrose export}}, 'numDel', max
% no. of deletions ie you can knockout a max of, say, 5 not more)

% Syntax - constrOpt = struct('rxnList', {{biomass}}, 'values', {{minimum
% accpetable biomass value ie lb}}, 'sense', 'G')
options = struct('targetRxn', char(model.rxns(1074)), 'numDel', 2);
constrOpt = struct('rxnList', {{'Biomass_Auto_2973'}},'values', 0.5*solWT.f,...
 'sense', 'G');


% This is to get the top "5" Knockouts (could be single knockouts, pairs, 
% triplets etc; change this number to how many ever differnet knockouts 
% you wanna get
no_of_knockouts = 5;

% These are to store all the knockouts that the code gives in each loop
% This way you can keep the code from coming up with the same knockouts
% each time
previousSolutions = cell(10, 1);

% 
contPreviousSolutions = 1;

% To find the knockouts
i = 1;

while i <= no_of_knockouts
    fprintf('...Performing optKnock analysis...\n')
    
    % If this is the first time we are performing OptKnock we do not have
    % to Worry about previous solutions
    if isempty(previousSolutions{1})
        optKnockSol = OptKnock(model, model.rxns, options, constrOpt);
        
    % If it is not then we have to keep in mind the reactions that were
    % already considered as potential knockouts
    else
        optKnockSol = OptKnock(model, model.rxns, options, constrOpt,...
            previousSolutions, 1);
    end
    
    ko_rxn = optKnockSol.rxnList;
    % Once you have performed OptKnock you need to save the values of
    % fluxes and the reaction that was deleted
    if ~isempty(ko_rxn)
        
        % Store the values
        previousSolutions{contPreviousSolutions} = ko_rxn;       
        contPreviousSolutions = contPreviousSolutions + 1;
        
        % Display all the reactions that have to be deleted
        fprintf('optKnock found a optKnock set of large %d composed by ',...
            length(ko_rxn))
        
        % To print the group of reactions to be deleted together to get the
        % above mentioned results
        for j = 1:length(ko_rxn)
            if j == 1
                fprintf('%s', ko_rxn{j})
            elseif j == length(ko_rxn)
                fprintf(' and %s', ko_rxn{j})
            else
                fprintf(', %s', ko_rxn{j})
            end
        end
    
        % This is the maximum possible yield disregarding every other reaction
        fprintf('\n')
        fprintf('The production of sucrose after optimization is %.2f \n',...
            optKnockSol.fluxes(1074))
        fprintf('The growth rate after optimization is %.2f \n', ...
            optKnockSol.fluxes(890))
    
    
        % OptKnock gives you all the possible knockouts to serve your purpose
        % But we will need something of a compormise between biomass and
        % sucrose production 
        fprintf('...Performing Coupling Analysis...\n')
        
     
        [type, maxGrowth, maxProd, minProd] = analyzeOptKnock(model, ko_rxn,...
            model.rxns(1074));
    
    
        fprintf('The solution is of type: %s\n', type);
        fprintf('The maximun growth rate given the optKnock set is %.2f\n',...
            maxGrowth)
        fprintf(['The maximun and minimun production of sucrose given the optKnock set is ' ...
            '%.2f and %.2f, respectively \n\n'], minProd, maxProd)
        
    else
        if i == 1
            fprintf('optKnock was not able to found an optKnock set\n')
        else
            fprintf('optKnock was not able to found additional optKnock sets\n')
        end
        break;
    end
    i = i + 1;
end
% ---------------------- END OF NO SALT STRESS ---------------------------
%}



%{
% --------------------------- SALT STRESS --------------------------------

% Conditions to keep the cells alive and for maximizing sucrose production
% reaction No. 890 is the biomass reaction and 1074 is sucrose export
model.lb(890) = 0.5 * solWT.v(890);
model.c(890)= 0;
model.c(1074) = 1;


% Some useful variables
val = zeros(2, length(SelRxn));
l = 0;
u = 0;



% For single knockouts 
for i = 1:length(SelRxn)

    % To Find the reaction number in the model
    for j = 1:length(model.rxns())
        if strcmp(char(model.rxns(j)), char(SelRxn(i)))
            index = j;
            break
        end
    end
    
    % These are to store is initial bounds so we can reverse the deletion
    % after getting results
    l = model.lb(index);
    u = model.ub(index);
    
    % Simulating Deletions
    model.lb(index) = 0;
    model.ub(index) = 0;
    
    % Find results for growth rates and sucrose production rates
    sol = optimizeCbModel(model);
    
    % To check if the solver was successful
    if sol.stat ~= 0
        val(1, i) = sol.v(890);
        val(2, i) = sol.v(1074);
    else 
        val(1, i) = 0;
        val(2, i) = 0;
    end
    
    % reversing the deletions
    model.lb(index) = l;
    model.ub(index) = u;
end


% Now we analyze the values for growth rate and sucrose production 
% all of which is stored in val
% This analysis will rank all the single knockouts based on the product of
% biomass and sucrose production
% Sort Method - Bubble Sort
for i = 1:length(SelRxn)
    
    for j = 1:length(SelRxn) - i
        
         cur = val(2, j) * val(1, j);
         nex = val(2, j+1) * val (1, j+1);
         
         
         if cur < nex
             
            temp = val(2, j);
            val(2, j) = val(2, j+1);
            val (2, j+1) = temp;
            
            temp = val(1, j);
            val(1, j) = val(1, j+1);
            val(1, j+1) = temp;
            
            temp = SelRxn(j);
            SelRxn(j) = SelRxn(j+1);
            SelRxn(j+1) = temp;
            
         end
   
    end
    
end

% Graphs
scatter3(1:length(SelRxn) , val(1, :), val(2, :), 5, 'filled')
xlabel('Reaction No.')
ylabel('Biomass')
zlabel('Sucrose')
% ------------------------- END OF SALT STRESS ---------------------------
%}
% ------------------  E N D  O F  O P T K N O C K-------------------------
%}



%{
% This variable keeps record of the target reaction and the limit on the
% number of knockouts 
% 
% Syntax - options = struct('targetRxn', {{butanol export}}, 'numDel', max
% no. of deletions ie you can knockout a max of, say, 5 not more)
% 
% Syntax - constrOpt = struct('rxnList', {{biomass}}, 'values', {{minimum
% acceptable biomass value ie lb}}, 'sense', 'G')
options = struct('targetRxn', char(model.rxns(1074)), 'numDel', 5);
constrOpt = struct('rxnList', {{'Biomass_Auto_2973'}},'values', 0.5*solWT.f,'sense', 'G');
selectedRxnList = {'PTAr';'PTA2'; 'ACALD'; 'ALCD2x'; 'ALCD19'; 'FRD2' ; 'FRD3';'LDH_D'; 'PGI'; 'PFK_3';'PFK_2';'PFK';'FBA';'FBA3'; 'TPI'; 'GAPD'; 'PGK'; 'PGM'; 'ENO'; 'PYK'; 'ACKr'; 'PFL'; 'OBTFL'; 'PDH'; 'ACONTb'; 'ACONTa' ; 'MICITDr'; 'ICDHyr'; 'AKGDH'; 'SUCOAS'; 'SUCDi'; 'FRD2'; 'FRD3' ; 'FUM'; 'DTARTD'; 'MDH'; 'CS'; 'PPC'}; 

threshold = 5;

previousSolutions = cell(10, 1);
contPreviousSolutions = 1;

% We will try to find 10 optKnock sets of a maximun length of 2
nIter = 1;
while nIter < threshold
    fprintf('...Performing optKnock analysis...')
    if isempty(previousSolutions{1})
        optKnockSol = OptKnock(model, model.rxns, options, constrOpt);
    else
        optKnockSol = OptKnock(model, model.rxns, options, constrOpt, previousSolutions);
    end
    % determine lactate production and growth rate after optimization
    sucFluxM1 = optKnockSol.fluxes(strcmp(model.rxns, char(model.rxns(1074))));
    growthRateM1 = optKnockSol.fluxes(strcmp(model.rxns,'Biomass_Auto_2973'));
    setM1 = optKnockSol.rxnList;
    if ~isempty(setM1)
        previousSolutions{contPreviousSolutions} = setM1;
        contPreviousSolutions = contPreviousSolutions + 1;
        %printing results
        fprintf('optKnock found a optKnock set of length %d composed by ', length(setM1));
        for j = 1:length(setM1)
            if j == 1
            fprintf('%s', setM1{j});
            elseif j == length(setM1)
            fprintf(' and %s', setM1{j});
            else
            fprintf(', %s', setM1{j});
            end
        end
        fprintf('\n');
        fprintf('The production of butanol after optimization is %.2f \n', sucFluxM1);
        fprintf('The growth rate after optimization is %.2f \n', growthRateM1);
        
        fprintf('...Performing coupling analysis...\n');
        
        [type, maxGrowth, maxProd, minProd] = analyzeOptKnock(model, setM1, char(model.rxns(1074)), 'Biomass_Auto_2973');
        fprintf('The solution is of type: %s\n', type);
        fprintf('The maximum growth rate given the optKnock set is %.2f\n', maxGrowth);
        fprintf(['The maximun and minimun production of butanol given the optKnock set is %.2f and %.2f, respectively \n\n'], minProd, maxProd);
        %singleProductionEnvelope(model, setM1, char(model.rxns(1074)), biomass, 'savePlot', 1, 'showPlot', 1, 'fileName', ['lact_ex2_' num2str(nIter)], 'outputFolder', 'OptKnockResults');
    else
    if nIter == 1
        fprintf('optKnock was not able to found an optKnock set\n');
    else
        fprintf('optKnock was not able to found additional optKnock sets\n');
    end
    break;
    end
    nIter = nIter + 1;
end
%}