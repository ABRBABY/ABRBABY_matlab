function [] = abrbaby_process_visual_comparison_by_group_age(eeglab_path, biosig_installer_path,indir,STD_number) 
% Load EEGLAB 
% addpath(genpath('/Users/anne-sophiedubarry/Documents/4_Software/eeglab2020_0'));
tmp = pwd ; 
cd(eeglab_path) ; 
% Open eeglab
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
run(biosig_installer_path) ; 
cd(tmp) ;

%Set group of electrodes to display for visualization
%Electrodes indices: 1='Fp1'; 2='Fp2';3='F4'; 4='Fz'; 5='F3'; 6='T7'; 7='C3'; 8='Cz'; 9='C4'; 
%10='T8'; 11='P4'; 12='Pz'; 13='P3'; 14='O1'; 15='Oz'; 16='O2'; 17='Lmon'; 
%18='Ref'; 19='Rmon'; 20='Left3'; 21='Right'
elec_to_disp_labels = {'F3','Fz','F4';'C3','Cz','C4'};
elec_indices = [5,4,3;7,8,9];

%Set 'STD_setfile' variable depending on STD_number (balanced or unbalanced number of STD trials)
if STD_number == 1
        STD_setfile = {('_STD1.set'),('_STD2.set')};
    elseif STD_number == 2
        STD_setfile = {('_all_STD1.set'),('_all_STD2.set')};
end
    
% Modify preferences in order to be able to load multiple datasets 
pop_editoptions( 'option_storedisk', 1);

% Reads all folders that are in indir 
subjects = filter_subjects_based_rejection(indir,0.2); 

%conditions = {'STD1', 'STD2', 'DEV1','DEV2'} ; 
%commands = {}; % initialize STUDY dataset list

grpA.suffix = {'_T3','_T6','_T8','_T10'};
grpB.suffix = {'_T18','_T24'};

%Get grand averages for each condition and subject
%(dimensions = subjects x channels x timepoints)
[grpA.DEV1_avg, grpA.DEV2_avg, grpA.STD1_avg, grpA.STD2_avg, timepoints, labels] = extract_averages_DEV_STD(indir,subjects(contains(subjects,grpA.suffix)),STD_setfile);
[grpB.DEV1_avg, grpB.DEV2_avg, grpB.STD1_avg, grpB.STD2_avg, timepoints, labels] = extract_averages_DEV_STD(indir,subjects(contains(subjects,grpB.suffix)),STD_setfile);

%Get standard error for each condition (dimensions = channels x timepoints)
[grpA.DEV1_se, grpA.DEV2_se,grpA.STD1_se, grpA.STD2_se] = extract_standard_error_DEV_STD(grpA.DEV1_avg, grpA.DEV2_avg, grpA.STD1_avg, grpA.STD2_avg)
[grpB.DEV1_se, grpB.DEV2_se,grpB.STD1_se, grpB.STD2_se] = extract_standard_error_DEV_STD(grpB.DEV1_avg, grpB.DEV2_avg, grpB.STD1_avg, grpB.STD2_avg)

%Mean activity through subjects (all electrodes)
grd_STD1_grpA = squeeze(mean(grpA.STD1_avg(:,:,:),1)) ; 
grd_STD2_grpA = squeeze(mean(grpA.STD2_avg(:,:,:),1)) ; 
grd_DEV1_grpA = squeeze(mean(grpA.DEV1_avg(:,:,:),1)) ; 
grd_DEV2_grpA = squeeze(mean(grpA.DEV2_avg(:,:,:),1)) ;  
grd_STD1_grpB = squeeze(mean(grpB.STD1_avg(:,:,:),1)) ; 
grd_STD2_grpB = squeeze(mean(grpB.STD2_avg(:,:,:),1)) ; 
grd_DEV1_grpB = squeeze(mean(grpB.DEV1_avg(:,:,:),1)) ; 
grd_DEV2_grpB = squeeze(mean(grpB.DEV2_avg(:,:,:),1)) ;

%Difference DEV1-STD1 for grpA
diff_DEV1_STD1_grpA = grd_DEV1_grpA - grd_STD1_grpA;

%Difference DEV1-STD1 for grp B
diff_DEV1_STD1_grpB = grd_DEV1_grpB - grd_STD1_grpB;

%Difference DEV2-STD2 for grpA
diff_DEV2_STD2_grpA = grd_DEV2_grpA - grd_STD2_grpA;

%Difference DEV2-STD2 for grpB
diff_DEV2_STD2_grpB = grd_DEV2_grpB - grd_STD2_grpB;

%Mean STD1+STD2 for grpA
mean_STD1_STD2_grpA = (grd_STD1_grpA + grd_STD2_grpA)/2;

%Mean STD1+STD2 for grpB
mean_STD1_STD2_grpB = (grd_STD1_grpB + grd_STD2_grpB)/2;


%% Visualisation

% %contains(labels,elec_to_disp_labels);  %[5,4,3;7,8,9];
% %elec_indices = find(contains(labels,elec_to_disp_labels)); %encours TODO
% %CODE TO WORK ON:
% %[sharedvals,idx] = intersect(labels, elec_to_disp_labels) ;
% %[sharedvals,other,idx] = intersect(elec_to_disp_labels,labels,'stable');
% %reshape(idx, size(elec_to_disp_labels)) ;
% 
% %try this :
% %------------------------------------------------------------------
% % elec= {'F3','Fz','F4';'C3','Cz','C4'};
% % labels = {'Fp1','Fp2','F4','Fz','F3','T7','C3','Cz','C4','T3'} ; 
% %  
% % [sharedvals,idxA,idxB] = intersect(labels, elec(:),'stable') ; 
% % [~,id1,id2] = intersect(elec(:), sharedvals,'stable') ; 
% % elec_indices = reshape(idxA(id2),size(elec)) ; % BINGO!
% %-------------------------------------------------------------------
% 
% elec_indices = zeros(size(elec_to_disp_labels,1),size(elec_to_disp_labels,2));
% for n = 1:size((elec_to_disp_labels),1)
%     [sharedvals,other,idx] = intersect(elec_to_disp_labels(n,:),labels,'stable');
%     elec_indices(n,:)=idx';
% end
% 
nb_young = size(grpA.DEV1_avg,1);
nb_old = size(grpB.DEV1_avg,1);

elec_indices_temp = zeros(size(elec_to_disp_labels,1),size(elec_to_disp_labels,2));
for n = 1:size((elec_to_disp_labels),1)
    [sharedvals,other,idx] = intersect(elec_to_disp_labels(n,:),labels,'stable');
    elec_indices_temp(n,:)=idx';
end

if STD_number == 1
    fig_name_STD = 'Response to standard /DA/ group comparison balanced number of STD';
    fig_name_DEV1 = 'Response to deviant /BA/ group comparison balanced number of STD';
    fig_name_DEV2 = 'Response to deviant /GA/ group comparison balanced number of STD';
    fig_name_MMN1 = 'MMN in response to deviant /BA/ group comparison balanced number of STD';
    fig_name_MMN2 = 'MMN in response to deviant /GA/ group comparison balanced number of STD';
    fig_name_cond1_grpA = 'MMN in response to deviant /BA/ 6-10mo group balanced number of STD';
    fig_name_cond1_grpB = 'MMN in response to deviant /BA/ 18-24mo group balanced number of STD';
    fig_name_cond2_grpA = 'MMN in response to deviant /GA/ 6-10mo group balanced number of STD';
    fig_name_cond2_grpB = 'MMN in response to deviant /GA/ 18-24mo group balanced number of STD';

elseif STD_number == 2
    fig_name_STD = 'Response to standard /DA/ group comparison unbalanced number of STD';
    fig_name_DEV1 = 'Response to deviant /BA/ group comparison unbalanced number of STD';
    fig_name_DEV2 = 'Response to deviant /GA/ group comparison unbalanced number of STD';
    fig_name_MMN1 = 'MMN in response to deviant /BA/ group comparison unbalanced number of STD';
    fig_name_MMN2 = 'MMN in response to deviant /GA/ group comparison unbalanced number of STD';
    fig_name_cond1_grpA = 'MMN in response to deviant /BA/ 6-10mo group unbalanced number of STD';
    fig_name_cond1_grpB = 'MMN in response to deviant /BA/ 18-24mo group unbalanced number of STD';
    fig_name_cond2_grpA = 'MMN in response to deviant /GA/ 6-10mo group unbalanced number of STD';
    fig_name_cond2_grpB = 'MMN in response to deviant /GA/ 18-24mo group unbalanced number of STD';
end

% figure; 
% plot([1,2,3],[1,2,3],'b'); hold on ; plot([3,2,1],[1,2,3],'r')
% leg1 = 'red line';
% leg2 = ['blue' newline 'line'];
% legend({leg1,leg2});

legend1_STD = ['Std /DA/ 6-10mo' newline ['n=', num2str(nb_young)]];
legend2_STD = ['Std /DA/ 18-24m' newline ['n=',num2str(nb_old)]];
plot_mean_STD = visual_comparison_several_electrodes_2values(elec_to_disp_labels,elec_indices,mean_STD1_STD2_grpA,mean_STD1_STD2_grpB,timepoints,fig_name_STD,legend1_STD,legend2_STD,STD_number);
%std_dev_patch(plot_mean_STD, grpA.STD1_se, grpB.STD1_se, color)

legend1_DEV1 = ['Dev /BA/ 6-10mo' newline ['n=', num2str(nb_young)]];
legend2_DEV1 = ['Dev /BA/ 18-24mo' newline ['n=', num2str(nb_old)]];
plot_mean_DEV1 = visual_comparison_several_electrodes_2values(elec_to_disp_labels,elec_indices,grd_DEV1_grpA,grd_DEV1_grpB,timepoints,fig_name_DEV1,legend1_DEV1,legend2_DEV1,STD_number);

legend1_DEV2 = ['Dev /GA/ 6-10mo'  newline ['n=', num2str(nb_young)]];
legend2_DEV2 = ['Dev /GA/ 18-24mo' newline ['n=', num2str(nb_young)]];
plot_mean_DEV2 = visual_comparison_several_electrodes_2values(elec_to_disp_labels,elec_indices,grd_DEV2_grpA,grd_DEV2_grpB,timepoints,fig_name_DEV2,legend1_DEV2,legend2_DEV2,STD_number);

legend1_MMN1 = ['MMN to /BA/ dev 6-10mo' newline ['n=',num2str(nb_young)]];
legend2_MMN1 = ['MMN to /BA/ dev 18-24mo' newline ['n=',num2str(nb_old)]];
plot_mean_MMN1 = visual_comparison_several_electrodes_2values(elec_to_disp_labels,elec_indices,diff_DEV1_STD1_grpA,diff_DEV1_STD1_grpB, timepoints,fig_name_MMN1,legend1_MMN1,legend2_MMN1,STD_number);

legend1_MMN2 = ['MMN to /GA/ dev 6-10mo' newline ['n=',num2str(nb_young)]];
legend2_MMN2 = ['MMN to /GA/ dev 18-24mo' newline ['n=',num2str(nb_old)]];
plot_mean_MMN2 = visual_comparison_several_electrodes_2values(elec_to_disp_labels,elec_indices,diff_DEV2_STD2_grpA,diff_DEV2_STD2_grpB, timepoints,fig_name_MMN2,legend1_MMN2,legend2_MMN2,STD_number);

legend1_cond1_grpA = ['Std /DA/ 6-10mo' newline ['n=',num2str(nb_young)]];
legend2_cond1_grpA = ['Dev /BA/ 6-10mo' newline ['n=',num2str(nb_young)]];
legend3_cond1_grpA = ['MMN to /BA/ dev 6-10mo, n=',num2str(nb_young)];
plot_mean_cond1_grpA = visual_comparison_several_electrodes_3values(1,elec_to_disp_labels,elec_indices,grd_STD1_grpA,grd_DEV1_grpA,diff_DEV1_STD1_grpA, timepoints,fig_name_cond1_grpA,legend1_cond1_grpA,legend2_cond1_grpA, legend3_cond1_grpA,STD_number);

legend1_cond1_grpB = ['Std /DA/ 18-24mo' newline ['n=',num2str(nb_old)]];
legend2_cond1_grpB = ['Dev /BA/ 18-24mo' newline ['n=',num2str(nb_old)]];
legend3_cond1_grpB = ['MMN to /BA/ dev 18-24mo' newline ['n=',num2str(nb_old)]];
plot_mean_cond1_grpB = visual_comparison_several_electrodes_3values(1,elec_to_disp_labels,elec_indices,grd_STD1_grpB,grd_DEV1_grpB,diff_DEV1_STD1_grpB, timepoints,fig_name_cond1_grpB,legend1_cond1_grpB,legend2_cond1_grpB, legend3_cond1_grpB,STD_number);

legend1_cond2_grpA = ['Std /DA/ 6-10mo' newline ['n=',num2str(nb_young)]];
legend2_cond2_grpA = ['Dev /GA/ 6-10mo' newline ['n=',num2str(nb_young)]];
legend3_cond2_grpA = ['MMN to /GA/ dev 6-10mo' newline ['n=',num2str(nb_young)]];
plot_mean_cond2_grpA = visual_comparison_several_electrodes_3values(2,elec_to_disp_labels,elec_indices,grd_STD2_grpA,grd_DEV2_grpA,diff_DEV2_STD2_grpA, timepoints,fig_name_cond2_grpA,legend1_cond2_grpA,legend2_cond2_grpA, legend3_cond2_grpA,STD_number);


legend1_cond2_grpB = ['Std /DA/ 18-24mo' newline ['n=',num2str(nb_old)]];
legend2_cond2_grpB = ['Dev /GA/ 18-24mo' newline ['n=',num2str(nb_old)]];
legend3_cond2_grpB = ['MMN to /GA/ dev 18-24mo' newline ['n=',num2str(nb_old)]];
plot_mean_cond2_grpB = visual_comparison_several_electrodes_3values(2,elec_to_disp_labels,elec_indices,grd_STD2_grpB,grd_DEV2_grpB,diff_DEV2_STD2_grpB, timepoints,fig_name_cond2_grpB,legend1_cond2_grpB,legend2_cond2_grpB, legend3_cond2_grpB,STD_number);



end
