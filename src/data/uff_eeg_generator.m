% Script to convert .mat EEG files into .unv format
%% DEAP Version
for subject = 1:32
    
    for trial= 1:40
        if subject <= 9
            load_name1=['s0',num2str(subject),'.mat']
        else
            load_name1=['s',num2str(subject),'.mat']
        end
        subject1=subject;
        s01=load(load_name1);
        basename=['S' num2str(subject1)];
        extension='.png';
        mkdir('export/uff', basename)
        
        
        Y1=s01.data(trial,(1:32),:);
        Y1=squeeze(Y1);
        t = linspace(0,63,length(Y1));
        fs=128;
        T=1/fs;
        
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
        
        filename=['Subject',num2str(subject),'_Trial',num2str(trial),'.uff']
        
        Info=writeuff(['export/uff/S',num2str(subject),'/',filename], UFF, 'replace');
    end
    
    
    
end
    
    
    
