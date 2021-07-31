%% Data extractor for pavlovia data
%
% Author: Filippo Cerpelloni
%
% IMPORTANTE! Nella cartella da cui si va ad estrarre (cd...) devono essere
% presenti solo i csv con i dati dei partecipanti, altrimenti genera un
% errore cercando di estrarre tutto

clear
% Cambiare path della cartella se i file vengono spostati
cd('C:\Users\filip\Dropbox\neuro\tirocinio post laurea\online\MPT_final\analisi_mpti\dati_mpti') 
datafile = dir('*.csv'); % tutti i csv salvati in una matrice
nData = length(datafile); % numero di csv = numero dati raccolti

% CAMBIARE PATH
cd('C:\Users\filip\Dropbox\neuro\tirocinio post laurea\online\MPT_final\analisi_mpti'); 

nPar = 1; 
for i = 1:nData 
    path = strcat(datafile(i).folder,'\',datafile(i).name); % concatena stringhe per percorso completo
    par_temp = mpti_import(path); % partecipante temporaneo
    % Se la lunghezza della tabella è 252, l'esperimento è stato svolto
    % completamente e il partecipante viene considerato, assegnandoli un
    % numero. Altrimenti, la tabella viene scartata.
    % (non ho trovato il modo di non estrarre affatto la tabella a priori)
    
    if size(par_temp,1) == 252
        eval(['par' num2str(nPar) ' = par_temp;']); % temporaneo diventa partecipante
        nPar = nPar+1;
    end
end

nPar = nPar -1;
%% Extraction of matrixes simplified for each participant 
% For each participant, from the original table that is conserved, a matrix
% with Yes / No / Total / Percentage is extracted to test for the curve
% fitting

for j = 1:nPar
    % eval evaluates the string as a command line. The expression would be: 
    % [sParj, dParj] = mpti_sort(parj);
    % Call to the function
    eval(['[sPar' num2str(j) ', dPar' num2str(j) '] = mpti_sort(par' num2str(j) ');']);
end

% salvataggio dati
save('confidence_an.mat');

%% Creation of matrixes for fitting analysis
% made to run on single sections

clear
load('confidence_an.mat'); % loads data extracted
clearvars datafile i j path nData

% Calcolo dei totali (somma delle risposte)
sTot = zeros(8,15);  
dTot = zeros(8,15);
for i = 1:nPar
    % Aggiunge al totale ogni partecipante (nPar da import .mat) per same 
    eval(['sPar = sPar' num2str(i) '(1:5,:); ']); % extracts the single results for each participants
    sTot(1:5,:) = sTot(1:5,:) + sPar;             % sum of the values into the total for Yes / No / 3
    % e diff condition
    eval(['dPar = dPar' num2str(i) '(1:5,:); ']); % extracts the single results for each participants
    dTot(1:5,:) = dTot(1:5,:) + dPar;             % sum of the values into the total for Yes / No / 3

end

% Calcolo dei totali (medie)
sTot(6,:) = sTot(1,:) + sTot(2,:) + sTot(3,:) + sTot(4,:) + sTot(5,:);
sTot(7,:) = sTot(1,:) ./ (sTot(1,:) + sTot(5,:));
sTot(8,:) = sTot(2,:) ./ (sTot(2,:) + sTot(4,:));

% per diff condition
dTot(6,:) = dTot(1,:) + dTot(2,:) + dTot(3,:)+ dTot(4,:) + dTot(5,:);
dTot(7,:) = dTot(1,:) ./ (dTot(1,:) + dTot(5,:));
dTot(8,:) = dTot(2,:) ./ (dTot(2,:) + dTot(4,:));

% tails for fitting curve (x axis)
SOA_L = [-350 -300 -250 -200 -150 -100 -50 0];
SOA_R = [0 50 100 150 200 250 300 350];

% same position, R and L tails (y axis)
sL = sTot(8,1:8); 
sR = sTot(8,8:15);

% different position, R and L tails (y axis)
dL = dTot(8,1:8);
dR = dTot(8,8:15);

save('confidence_an_input.mat');

%% Curve fitting script

clear
load('confidence_an_input.mat');

% createFits: last parameter (logical) indicates whether to visualize or
% not the graphs (in the single participants, it's 4 graph for each one, it
% could be a lot)

% GOF indexes:              
%      Position  Tail
% sL = same      left 
% sR = same      right
% dL = diff      left
% dR = diff      right

% General curve fitting (with plot)
%[results, gof] = mpti_createFits(SOA_L, sL, SOA_R, sR, dL, dR, true);

% Final table initialized
toWrite = table; 

% Single participant fittings
% - add participant details
% - get fitting data for same, diff, sure, maybe, every SOA
% - put everything in the final table

for j = 1:nPar
    % Get temp values, eval is not a great function
    % e.g. sl_temp = sPar1(5,1:8);
    % Postision-Certainty-Tail: e.g. smr = same-maybe-right
    eval(['ssr_temp = sPar' num2str(j) '(7,1:8);']);  eval(['ssl_temp = sPar' num2str(j) '(7,8:15);']);   
    eval(['smr_temp = sPar' num2str(j) '(8,1:8);']);  eval(['sml_temp = sPar' num2str(j) '(8,8:15);']);
    eval(['dsr_temp = sPar' num2str(j) '(7,1:8);']);  eval(['dsl_temp = sPar' num2str(j) '(7,8:15);']);   
    eval(['dmr_temp = sPar' num2str(j) '(8,1:8);']);  eval(['dml_temp = sPar' num2str(j) '(8,8:15);']);
        
    % Curve fitting salvati in variabili resJ e gofJ
    eval(['[resS' num2str(j) ', gofS' num2str(j) '] = mpti_createFits(SOA_L, ssl_temp, SOA_R, ssr_temp, dsl_temp, dsr_temp, false);']);
    eval(['[resM' num2str(j) ', gofM' num2str(j) '] = mpti_createFits(SOA_L, sml_temp, SOA_R, smr_temp, dml_temp, dmr_temp, false);']);
    
    % Temp vars to ease next steps
    eval(['rs = resS' num2str(j) ';']); % r = resJ;
    eval(['gs = gofS' num2str(j) ';']); % g = gofJ;
    eval(['rm = resM' num2str(j) ';']); 
    eval(['gm = gofM' num2str(j) ';']); 
    eval(['p = par' num2str(j) '(1,:);']); % prima riga della tabella partecipante
    eval(['s = sPar' num2str(j) ';']); % s = sParJ; Per inserire i rate grezzi
    eval(['d = dPar' num2str(j) ';']); % d = dParJ;

    % temp table, to write less letters 
    t = table; 
    
    % Add data from participants
    t.Iniziali = string(p.InizialinomeecognomeesMarioRossiMR);
    t.Data_nascita = p.DatadinascitaGGMMAAAA;
    t.Genere = p.GenereMF;
    t.Mano_dominante = p.ManodominanteDS;
    t.Scolarita = p.Anniscolaritdiploma13laureatriennale16etc;
    t.Data = p.date;
    t.PsychoPy_version = p.psychopyVersion;
    t.OS = p.OS;
    t.Frame_rate = p.frameRate;
    
    % Add raw ratings data, organized per congurency, certainty, SOA 
    for soa = 1:15
        eval(['t.Rate_Same_Sure_SOA' num2str(soa) '  = s(7,' num2str(soa) ');']);
    end
    for soa = 1:15
        eval(['t.Rate_Same_Maybe_SOA' num2str(soa) '  = s(8,' num2str(soa) ');']);
    end
    for doa = 1:15
        eval(['t.Rate_Diff_Sure_SOA' num2str(doa) '  = d(7,' num2str(doa) ');']);
    end
    for doa = 1:15
        eval(['t.Rate_Diff_Maybe_SOA' num2str(doa) '  = d(8,' num2str(doa) ');']);
    end
      
    % Add fittings data: slope, thresholds 
    % (50 and 75, both value and absolute), adjustedR 
    % Organized per side, condition (Audio or Visual Leading) 
    
    % Same position, Certainty high, Audio Leading
    t.Same_Sure_AL_slope = rs{1,1}.b;       t.Same_Sure_AL_soglia50 = rs{1,1}.t;   
    t.Same_Sure_AL_s50ABS = abs(rs{1,1}.t); t.Same_Sure_AL_adjRsquared = gs(1).adjrsquare;
    
    % Same position, Certainty low, Audio Leading
    t.Same_Maybe_AL_slope = rm{1,1}.b;       t.Same_Maybe_AL_soglia50 = rm{1,1}.t;   
    t.Same_Maybe_AL_s50ABS = abs(rm{1,1}.t); t.Same_Maybe_AL_adjRsquared = gm(1).adjrsquare;
    
    % Same position, Certainty high, Visual Leading
    t.Same_Sure_VL_slope = rs{2,1}.b;       t.Same_Sure_VL_soglia50 = rs{2,1}.t;   
    t.Same_Sure_VL_s50ABS = abs(rs{2,1}.t); t.Same_Sure_VL_adjRsquared = gs(2).adjrsquare;
    
    % Same position, Certainty low, Visual Leading
    t.Same_Maybe_VL_slope = rm{2,1}.b;       t.Same_Maybe_VL_soglia50 = rm{2,1}.t;   
    t.Same_Maybe_VL_s50ABS = abs(rm{2,1}.t); t.Same_Maybe_VL_adjRsquared = gm(2).adjrsquare;
    
    % TBW for each threshold and certainty
    t.TBW_SAME_SURE_50 = abs(rs{1,1}.t)+abs(rs{2,1}.t); 
    t.TBW_DIFF_MAYBE_50 = abs(rm{1,1}.t)+abs(rm{2,1}.t);
    
    % Diff position, Certainty high, Audio Leading
    t.Diff_Sure_AL_slope = rs{3,1}.b;       t.Diff_Sure_AL_soglia50 = rs{3,1}.t;   
    t.Diff_Sure_AL_s50ABS = abs(rs{3,1}.t); t.Diff_Sure_AL_adjRsquared = gs(3).adjrsquare;
    
    % Diff position, Certainty low, Audio Leading
    t.Diff_Maybe_AL_slope = rm{3,1}.b;       t.Diff_Maybe_AL_soglia50 = rm{3,1}.t;   
    t.Diff_Maybe_AL_s50ABS = abs(rm{3,1}.t); t.Diff_Maybe_AL_adjRsquared = gm(3).adjrsquare;
       
    % Diff position, Certainty high, Visual Leading
    t.Diff_Sure_VL_slope = rs{4,1}.b;       t.Diff_Sure_VL_soglia50 = rs{4,1}.t;   
    t.Diff_Sure_VL_s50ABS = abs(rs{4,1}.t); t.Diff_Sure_VL_adjRsquared = gs(4).adjrsquare;
    
    % Diff position, Certainty low, Visual Leading
    t.Diff_Maybe_VL_slope = rm{4,1}.b;       t.Diff_Maybe_VL_soglia50 = rm{4,1}.t;   
    t.Diff_Maybe_VL_s50ABS = abs(rm{4,1}.t); t.Diff_Maybe_VL_adjRsquared = gm(4).adjrsquare;
    
    % TBW for each threshold and certainty
    t.TBW_DIFF_SURE_50 = abs(rs{3,1}.t)+abs(rs{4,1}.t); 
    t.TBW_DIFF_MAYBE_50 = abs(rm{3,1}.t)+abs(rm{4,1}.t);
    
    % Values across position
    % Join same and diff
    eval(['al = sPar' num2str(j) '(1:5,1:8) + dPar' num2str(j) '(1:5,1:8);']);
    al(6,:) = al(1,:) + al(2,:) + al(4,:) + al(5,:);
    al(7,:) = al(1,:) ./ (al(1,:) + al(5,:));
    al(8,:) = al(2,:) ./ (al(2,:) + al(4,:));
 
    eval(['vl = sPar' num2str(j) '(1:5,8:15) + dPar' num2str(j) '(1:5,8:15);']);
    vl(6,:) = vl(1,:) + vl(2,:) + vl(4,:) + vl(5,:);
    vl(7,:) = vl(1,:) ./ (vl(1,:) + vl(5,:));
    vl(8,:) = vl(2,:) ./ (vl(2,:) + vl(4,:));
    
    % createFits
    [rAll,gAll] = mpti_createFits(SOA_L, al(7,:), SOA_R, vl(7,:), al(8,:), vl(8,:), false);
        
    % Audio Leading
    t.AL_medio_Sure_slope = rAll{1,1}.b;          t.AL_medio_Sure_soglia50 = rAll{1,1}.t;
    t.AL_medio_Sure_s50ABS = abs(rAll{1,1}.t);    t.AL_medio_Sure_adjrsquared = gAll(1).adjrsquare;
    
    t.AL_medio_Maybe_slope = rAll{3,1}.b;         t.AL_medio_Maybe_soglia50 = rAll{3,1}.t;
    t.AL_medio_Maybe_s50ABS = abs(rAll{3,1}.t);   t.AL_medio_Maybe_adjrsquared = gAll(3).adjrsquare;   

    % Visual leading
    t.VL_medio_Sure_slope = rAll{2,1}.b;          t.VL_medio_Sure_soglia50 = rAll{2,1}.t;
    t.VL_medio_Sure_s50ABS = abs(rAll{2,1}.t);    t.VL_medio_Sure_adjrsquared = gAll(2).adjrsquare;
    
    t.VL_medio_Maybe_slope = rAll{4,1}.b;         t.VL_medio_Maybe_soglia50 = rAll{4,1}.t;
    t.VL_medio_Maybe_s50ABS = abs(rAll{4,1}.t);   t.VL_medio_Maybe_adjRsquared = gAll(4).adjrsquare;

    % TBW for each certainty and level
    t.TBW_TOT_SURE_50 = abs(rAll{1,1}.t) + abs(rAll{2,1}.t);  
    t.TBW_TOT_MAYBE_50 = abs(rAll{3,1}.t) + abs(rAll{4,1}.t);  
 
    % Last addition: out of the total, how many answers were sure and maybe
    eval(['sumS = sum(sPar' num2str(j) '([1 5],:) + dPar' num2str(j) '([1 5],:),''All'');']);
    eval(['sumM = sum(sPar' num2str(j) '([2 4],:) + dPar' num2str(j) '([2 4],:),''All'');']);
    sumT = sumM + sumS;
    
    t.Rate_Sure = sumS / sumT;
    t.Rate_Maybe = sumM / sumT;

    % Finally, add this row to the final table
    toWrite = [toWrite; t];
    
    % Data format checks
    % Whether participants put the correct name format, otherwise correct it 
    % e.g. 'Filippo Cerpelloni' -> FC
    a = strsplit(string(toWrite.Iniziali(j)),' ');
    if length(a) > 1  
        toWrite.Iniziali(j) = string(upper(a{1}(1))) + string(upper(a{2}(1))); % due iniziali maiuscole
    end
    
    % Check the birth date format, must be DDMMYYYY
    b = strsplit(string(toWrite.Data_nascita(j)),{'/','.','-'});
    if length(b) > 1  
        toWrite.Data_nascita(j) = string(b{1}) + string(b{2}) + string(b{3}); % concatena giorni mese anno
    end
    
    % Check that handedness and gender and in uppercase
    if isstrprop(string(toWrite.Mano_dominante(j)), 'lower') % mano dominante
        toWrite.Mano_dominante(j) = upper(string(toWrite.Mano_dominante(j)));
    end
    if isstrprop(string(toWrite.Genere(j)), 'lower') % genere
        toWrite.Genere(j) = upper(string(toWrite.Genere(j)));
    end
end

% !!! CAMBIAMENTI MANUALI !!!
% Eventualmente implementabili nel codice più avanti
toWrite = sortrows(toWrite,{'Iniziali','Data'},'ascend');
toWrite(end,:) = [];
toWrite.Data_nascita(14) = '28041997'; % Lettera extra nella data
toWrite.Data_nascita(15) = '13111997'; % Novembre scritto in parola
toWrite.Scolarita(35) = 19; % JP: Aggiunta di nome del titolo: master / diploma
toWrite.Scolarita(15) = 13;
toWrite.Scolarita(22) = 13; % FB 97

% salva file mat con tutti i dataset per e con le analisi
% prima, rimosse variabili temporanee
clearvars dL dl_temp dR dr_temp g i j a b p r sL sl_temp sR sr_temp t
save('mpti_fit_results.mat');

% Wirte xlsx and csv files  
fn_xlsx = "mpti_confidence_table.xlsx";  writetable(toWrite,fn_xlsx);
fn_csv = "mpti_risultati_n" + nPar + ".csv";    writetable(toWrite,fn_csv,'Delimiter',',');



