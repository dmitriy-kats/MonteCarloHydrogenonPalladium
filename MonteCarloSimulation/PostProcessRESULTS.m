clc; clear; close all;

FOLDERNUMBER=8;

for j=1:1
    mu_atfixedT=zeros(1,51);
    Theta_atfixedT=zeros(1,51);
    
    for i=0:1:200
    filename1=['./RESULTS' num2str(FOLDERNUMBER) '/ThetaValues' num2str(i+201*(j-1)) '.out']
    fileID1 = fopen(filename1);
    Theta = fread(fileID1, 'double');
    fclose('all');
    
    filename2=['./RESULTS' num2str(FOLDERNUMBER) '/inputs' num2str(i+201*(j-1)) '.txt'];
    fileID2 = fopen(filename2);
    cellData=textscan(fileID2,'%s %f',2,'Delimiter',',');
    fclose('all');
    
    filename3=['./RESULTS' num2str(FOLDERNUMBER) '/Cv' num2str(i+201*(j-1)) '.out'];
    fileID3 = fopen(filename3);
    Cv = fread(fileID3, 'double');
    fclose('all');
    
    
    muandTempInputs=cell2mat(cellData(2));
    %disp([i muandTempInputs(1) muandTempInputs(2) mean(Theta(900:end)) std(Theta)*1e16]);
    mu_atfixedT(i+1)=muandTempInputs(1);
    Theta_atfixedT(i+1)=mean(Theta);
    Cv_atfixedT(i+1)=mean(Cv);

    end
    
    figure
    plot(mu_atfixedT, Cv_atfixedT, '.', 'MarkerSize',8);
    xlabel('$\mu$', 'Interpreter', 'Latex', 'FontSize', 16);
    ylabel('Cv', 'Interpreter', 'Latex', 'FontSize', 16);
    titletext=['T:' num2str(muandTempInputs(2))];
    title(titletext, 'Interpreter', 'Latex', 'FontSize', 16);
    %axis([3 6 0.5 1.0])
    %pause 
    %close all
    
end


aa=[mu_atfixedT; Cv_atfixedT; Theta_atfixedT];