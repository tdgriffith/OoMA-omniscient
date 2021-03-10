% Script to convert .mat EEG files into .unv format
%% SPIS Version
for trial = 1:2
    if trial == 1
        test='EC';
    elseif trial == 2
        test = 'EO';
    end
    
    for subject= 2:11
        if subject <= 9
            load_name1=['S0',num2str(subject),'_restingPre_',test,'.mat']
        else
            load_name1=['S',num2str(subject),'_restingPre_',test,'.mat']
        end
        subject1=subject;
        s01=load(load_name1);
        basename=['S' num2str(subject1)];
        extension='.png';
        mkdir('export/uff_SPIS', basename)
        fs=256;
        T=1/fs;
        
        Y1=s01.dataRest((1:64),:);
        Y1=squeeze(Y1);
        Y1_filt=bandpass(detrend(Y1'),[3 50],fs);
        Y1_filt=Y1_filt';
        t = linspace(0,60*2.5,length(Y1_filt));
        
        
        data1 = iddata(Y1_filt',[],T);
        data1.UserData=[];
        
        
        for i=1:size(Y1_filt,1)
            data1.OutputName{i}=[num2str(i),'+Y'];
            data1.UserData.OutputNodeID{i}=i;
            data1.OutputUnits{i}='m/s^2';
        end
        
        UFF=IdData2UNI(data1);
        for i = 1:length(UFF)
            UFF{1,i}.x=transpose(UFF{1,i}.x);
            UFF{1,i}.measData=transpose(UFF{1,i}.measData);
            UFF{1,i}.abscUnitsLabel='s';
            UFF{1,i}.ordDataChar=12;
            UFF{1,i}.ordinateAxisLabel='Acceleration        ';
            UFF{1,i}.rspDir=3;
            UFF{1,i}.binary=0;
            UFF{1,i}.d2='NONE';
        end
        
        filename=['Subject',num2str(subject),'_Trial',test,'.uff']
        
        Info=writeuff(['export/uff_SPIS/S',num2str(subject),'/',filename], UFF, 'replace');
    
    end
    
    
    
end
    
    
    

