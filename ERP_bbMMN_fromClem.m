    eeglab

    clear all
    clear all


    % s=0;   1  2 3 4 5  6     8  9  10   11 12      14 15 16    19  20 23  24 25 26   27: 29 30 31 32 35:38
    suj=0; %1 2 3 4 5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38
     for w= [3 4 8 9 10 13 14 15 16 28 29  32 33  40  44 46 53  79 84 86]
        %        if ww~=2 && ww~=4 && ww~=5 && ww ~= 6 && ww ~= 8 && ww ~= 11 && ww ~= 23 && ww ~= 26

        %             s=s+1
        suj=suj+1
       if w<10
            EEG = pop_loadbv(['C:\Users\clem\Desktop\MMN_T0\'],['bb0' num2str(w) '_T0_MMN.vhdr']);
        elseif w>=10
            EEG = pop_loadbv(['C:\Users\clem\Desktop\MMN_T0\'],['bb' num2str(w) '_T0_MMN.vhdr']);
        end
        
        
        %%%%%%% TF parameters %%%%%%

        
            samplingrate=256;
            freq=1:40;
            SCALE=samplingrate*1.5./freq;

            %%%%%%%%%%%%%%%%%%FIltering%%%%%%%%%%%%%%%

                                     EEG  = pop_basicfilter( EEG,  1:16 , 'Cutoff', [1 20], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  4 ); 
                      EEG  = pop_basicfilter( EEG,  1:16 , 'Cutoff',  50, 'Design', 'notch', 'Filter', 'PMnotch', 'Order',  180, 'RemoveDC', 'on' ); 
%   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%Definition des matrices de données%%%%%%%%%%

            STDD=zeros(16,151);kA2l1=0;

           %%%%%
            CONS=zeros(16,151);kA3l1=0;

           %%%%%
            FREQ=zeros(16,151);kA4l1=0;

           %%%%%
            VOY=zeros(16,151);kA5l1=0;


            %%%%%%%%%%%%%%%%%  extraction des Standards et Deviants %%%%%%%%%%%%%%

            numtr=0;
            rej=0;
            numtr2=0;
            rej2=0;
            numtr3=0;
            rej3=0;
            numtr4=0;
            rej4=0;

            for i=3:size(EEG.event,2)-4

                 if strcmp(EEG.event(i).type,'65283')
                  AUX=EEG.data(:,EEG.event(i).latency+17:EEG.event(i).latency+167);

               else
                  AUX=EEG.data(:,EEG.event(i).latency-12:EEG.event(i).latency+138);
               end 

               AUX=AUX-repmat(mean(AUX(:,1:25),2),1,151); %Treure línia base

                %%%%
    %             TF=EEG.data(:,(EEG.event(i).latency-500):(EEG.event(i).latency+500));
    %           %  TF(17:20,:)=[];
    %            % TF=TF-repmat(mean(TF(17:18,:),1),18,1);  %%Treiem mastoides
    %             TF=TF-repmat(mean(TF(:,487:500),2),1,1001); %Treure línia base
    %             %%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%STANDARDS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if strcmp(EEG.event(i).type,'65281') 
                    numtr=numtr+1;

                    if numtr<=900

                        rej=rej+1;

                        if max(max(abs(AUX)))<120

                            kA2l1=kA2l1+1;
                            STDD=STDD + AUX;

                            REJ(rej)=1;

                        elseif max(max(abs(AUX)))>=120
                            REJ(rej)=0;

                        end
                    end

                    %%%%%%%%%%%%%DEVIANT CONSONNE%%%%%%%%%%%%%%%%

                elseif strcmp(EEG.event(i).type,'65283')
                    numtr2=numtr2+1;

                    if numtr2<=900

                        rej2=rej2+1;

                        if max(max(abs(AUX)))<120

                            kA3l1=kA3l1+1;
                          CONS=CONS + AUX;

                            REJ2(rej2)=1;


                      elseif max(max(abs(AUX)))>=120
                            REJ2(rej2)=0;

                        end
                    end


                     %%%%%%%%%%%%%DEVIANT FREQUENCE%%%%%%%%%%%%%%%%


                elseif strcmp(EEG.event(i).type,'65285')
                    numtr3=numtr3+1;

                    if numtr3<=900

                        rej3=rej3+1;

                        if max(max(abs(AUX)))<120

                            kA4l1=kA4l1+1;
                          FREQ=FREQ + AUX;

                            REJ3(rej3)=1;

                         elseif max(max(abs(AUX)))>=120
                            REJ3(rej3)=0;



                        end
                    end    

                %%%%%%%%%%%%%DEVIANT VOYELLE%%%%%%%%%%%%%%%%%%%%%%

                   elseif strcmp(EEG.event(i).type,'65287')
                    numtr4=numtr4+1;

                    if numtr4<=900

                        rej4=rej4+1;

                        if max(max(abs(AUX)))<120

                            kA5l1=kA5l1+1;
                          VOY=VOY + AUX;

                            REJ4(rej4)=1;


                      elseif max(max(abs(AUX)))>=120
                            REJ4(rej4)=0;


                        end
                    end    

                end
             end


            ALL_STDD(:,:,suj)=STDD./kA2l1; A2KL1(suj)=kA2l1;
            A2KREJ(suj)=mean(REJ);


            ALL_CONS(:,:,suj)=CONS./kA3l1; A3KL1(suj)=kA3l1;
            A3KREJ(suj)=mean(REJ2);

            ALL_FREQ(:,:,suj)=FREQ./kA4l1; A4KL1(suj)=kA4l1;
            A4KREJ(suj)=mean(REJ3);

            ALL_VOY(:,:,suj)=VOY./kA5l1; A5KL1(suj)=kA5l1;
            A5KREJ(suj)=mean(REJ4);

            %%%%



        end



%%%%%% TF %%%%%%%

MGTFSTDD_FT_LP10_T0=mean(ALL_STDD(:,:,:),3);
MGTFCONS_FT_LP10_T0=mean(ALL_CONS(:,:,:),3);
MGTFFREQ_FT_LP10_T0=mean(ALL_FREQ(:,:,:),3);
MGTFVOY_FT_LP10_T0=mean(ALL_VOY(:,:,:),3);

 
MMNCons_FT_LP10_T0=MGTFCONS_FT_LP10_T0-MGTFSTDD_FT_LP10_T0;
MMNFreq_FT_LP10_T0=MGTFFREQ_FT_LP10_T0-MGTFSTDD_FT_LP10_T0;
MMNVoy_FT_LP10_T0=MGTFVOY_FT_LP10_T0-MGTFSTDD_FT_LP10_T0;






