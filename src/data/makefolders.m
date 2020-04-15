for i=1:40
    basename=['Trials/T' num2str(i)];
    mkdir(basename)
end

%%

for i=1:32
    for ii=1:40
        filename=['S' num2str(i) '/' 'S' num2str(i) 'E' num2str(ii) '.png']
        copyfile(filename, ['Trials/T' num2str(ii)])
    end
end

