% Script to convert .mat EEG files into .unv format
%% SPIS Version
for trial = 1:4    
    for subject= 1:1
        if subject <= 9
            load_name1=['S',num2str(subject),'T0',num2str(trial),'.h5']
        else
            load_name1=['S',num2str(subject),'T0',num2str(trial),'.h5']
        end
        subject1=subject;
        s01 = h5read(load_name1,['/df/block0_values']);
        basename=['S' num2str(subject1)];
        extension='.png';
        mkdir('export/uff_seizure', basename)
        fs=256;
        T=1/fs;
        
        Y1=s01;
        Y1=squeeze(Y1);        
        
        data1 = iddata(Y1',[],T);
        data1.UserData=[];
        
        
        for i=1:size(Y1,1)
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
        
        filename=['Subject',num2str(subject),'_Trial',num2str(trial),'seizureMIT.uff'];
        
        Info=writeuff(['export/uff_seizure/S',num2str(subject),'/',filename], UFF, 'replace');
    
    end
    
    
    
end
    
    
    

