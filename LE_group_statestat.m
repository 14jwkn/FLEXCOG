%{
For the specified k, find the frequency metrics for 
each run separately, then average across runs. This includes occurrence, 
dwell time, transition number, and transition probability.
Adapts code from: 
https://github.com/trendscenter/gift
Output:
avg_substat.csv Frequency metrics averaged across runs for each subject. 
%}

%Define command line arguments.
function [] = LE_group_statestat(k)
disp(append('Doing: ',k))

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set up I/O.
subfile = 'r_full_submain.txt'; 
subgroup = 'full';
outpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                subgroup,'/',k,'/');

%If the output folder is not created, create it.
if not(isfolder(outpath))
    mkdir(outpath)
end

%Read in subjects and set parameters.
subjects = textread(subfile,'%s','delimiter',',');
subjects = string(subjects); 
nsubj = size(subjects,1);
k = str2num(k);

%Read in the clustering.
infile = append(outpath,'uni_subcent.h5');
inkey = append('/',subgroup,'/table');
fullclust = double(h5read(infile,inkey).values);

%Extract the unique values.
klabs = unique(fullclust)';

%Make stat labels.
stat_labs = {};
i = 1;
for j = klabs
    stat_labs{i} = append('occur_',num2str(j));
    i = i + 1;
end
for j = klabs
    stat_labs{i} = append('dwell_',num2str(j));
    i = i + 1;
end
stat_labs{i} = append('numtrans');
i = i + 1;
for j = klabs
    for l = klabs
        stat_labs{i} = append('sptrans_',num2str(j),'-',num2str(l));
        i = i + 1;
    end    
end

%Refine stat labels.
stat_labs = string(stat_labs);
stat_idx = [subjects].';

%Create run labels.
run_labs = {'r1_lr','r1_rl','r2_lr','r2_rl'};

%Set lengths for subject and run.
subwin = length(fullclust)/nsubj;
runwin = subwin/4;
Nwin = runwin;

%Produce index labels for the k labels.
indclust = zeros(length(fullclust),1);
for i = 1:k
    indclust(fullclust==klabs(i))=i;
end   

%For each run.
for runin = 1:4

    %Create empty stat matrix.
    statmat = zeros((nsubj),size(stat_labs,2));
    statmat = array2table(statmat,'RowNames',stat_idx,'VariableNames',stat_labs);

    %For each subject.
    for subin = 1:nsubj

        %Extract subject clustering.
        substart = uint64((subwin*(subin-1)) + 1);
        subend = uint64(subwin*subin);
        subfull = indclust(substart:subend);
        
        %Extract current run.
        runstart = uint64((runwin*(runin-1)) + 1);
        runend = uint64(runwin*runin);
        a = subfull(runstart:runend);

        %Fraction of time spent in each state.
        F = zeros(1,k);
        for jj = 1:k
            F(jj) = (sum(a == jj))/Nwin;
        end

        %Number of Transitions
        NT = sum(abs(diff(a)) > 0);

        %Mean dwell time in each state
        MDT = zeros(1,k);
        for jj = 1:k
            start_t = find(diff(a==jj) == 1);
            end_t = find(diff(a==jj) == -1);
            if a(1)==jj
                start_t = [0; start_t];
            end
            if a(end) == jj
                end_t = [end_t; Nwin];
            end
            MDT(jj) = mean(end_t-start_t);
            if isempty(end_t) & isempty(start_t)
                MDT(jj) = 0;
            end
        end

        %Full Transition Matrix
        TM = zeros(k,k);
        for t = 2:Nwin
            TM(a(t-1),a(t)) =  TM(a(t-1),a(t)) + 1;
        end
        for jj = 1:k
            if sum(TM(jj,:)>0)
                TM(jj,:) = TM(jj,:)/sum(a(1:Nwin-1) == jj);
            else
                TM(jj,jj) = 1;
            end
        end 

        %Input into stat matrix.
        sub = subjects(subin);
        startin = 1;
        endin = k;
        statmat(sub,startin:endin) = num2cell(F);
        startin = endin + 1;
        endin = startin + k - 1;
        statmat(sub,startin:endin) = num2cell(MDT);
        endin = endin + 1;
        statmat(sub,endin) = num2cell(NT);
        startin = endin + 1;
        endin = startin + (k*k) - 1;
        TM_vec = {};
        i = 1;
        for j = 1:k
            for l = 1:k
                TM_vec{i} = TM(j,l);
                i = i + 1;
            end    
        end
        statmat(sub,startin:endin) = TM_vec;
    end
    
    %Save.
    rlab = string(run_labs(runin));
    statfile = append(outpath,rlab,'_substat.csv');
    writetable(statmat,statfile,'WriteRowNames',1)
    disp('Saved.')
end

%Read in the run files and average.
statmat = zeros((nsubj),size(stat_labs,2));
for runin = 1:4
    
    %Read matrix.
    rlab = string(run_labs(runin));
    infile = append(outpath,rlab,'_substat.csv');
    runmat = readmatrix(infile);
    runmat(:,1) = [];
    
    %Add to empty matrix.
    statmat = statmat + runmat;
end
statmat = statmat / 4;

%Save average.
statmat = array2table(statmat,'RowNames',stat_idx,'VariableNames',stat_labs);
statfile = append(outpath,'avg_substat.csv');
writetable(statmat,statfile,'WriteRowNames',1)
disp('All saved.')
end
