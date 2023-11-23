%{
Produce sample dFC(t) and LE(t) to visualize for the methods description.
Output:
dFC.csv Contains the sample dFC(t).
LE.csv Contains the sample LE(t).
%}

%Define command line arguments.
function [] = methods(subject)

%Set up I/O.
infile = append('../outputs/r_meants/',subject,...
                '/demean_rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv');

%Set up parameters.
tc = readmatrix(infile)';
[all_ts, N_areas] = size(tc);

%Get the BOLD phase using the Hilbert transform.
BOLD = tc.';
Phase_BOLD=zeros(N_areas,all_ts);
for seed=1:N_areas
    BOLD(seed,:) = BOLD(seed,:)-mean(BOLD(seed,:)); %Demean.
    Phase_BOLD(seed,:) = angle(hilbert(BOLD(seed,:))); %Hilbert transform.
end
clear tc

%Select time points.
t_list = [11 22 33];
for t=t_list
    
    %Calculate the instantaneous phase coherence, iFC.
    iFC=zeros(N_areas);
    for n=1:N_areas
        for p=1:N_areas
            iFC(n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
        end
    end
    dFC = iFC;
    
    %Calculate the leading eigenvector.
    [V1,~]=eigs(iFC,1);
    if mean(V1>0)>.5
        V1=-V1;
    elseif mean(V1>0)==.5 && sum(V1(V1>0))>-sum(V1(V1<0))
        V1=-V1;
    end
    Leading_Eig=V1;
    
    %Save.
    outfile1 = append('../outputs/outcollect/dFC_',subject,'_',num2str(t),'.csv');
    outfile2 = append('../outputs/outcollect/LE_',subject,'_',num2str(t),'.csv');
    writematrix(dFC,outfile1)
    writematrix(Leading_Eig,outfile2)
end
end
