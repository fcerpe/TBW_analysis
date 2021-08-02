%% single plots - 1 participant: old plot v. split sure-maybe for both same and different sides

SOA=[-350 -300 -250 -200 -150 -100 -50 0 50 100 150 200 250 300 350];

for i = 1:56
    f = figure;
    f.Position = [50 50 900 700];
    % Calc old ratings 
    eval(['sPar' num2str(i) '(9,:) = (sPar' num2str(i) '(1,:)+sPar' num2str(i) '(2,:)) ./ sPar' num2str(i) '(6,:)']);
    eval(['dPar' num2str(i) '(9,:) = (dPar' num2str(i) '(1,:)+dPar' num2str(i) '(2,:)) ./ dPar' num2str(i) '(6,:)']);

    for sp = 1:6
        subplot(2,3,sp);
        switch sp 
            case 1
                eval(['plot(SOA, sPar' num2str(i) '(9,:), ''b-*'')']); % old rating: same position
                title("Par " + i + " - same side: joined ratings");
            case 2
                eval(['plot(SOA, sPar' num2str(i) '(7,:), ''b-*'')']); % old rating: same position
                title("Par " + i + " - same: certain ratings");
            case 3
                eval(['plot(SOA, sPar' num2str(i) '(8,:), ''b-*'')']); % old rating: same position
                title("Par " + i + " - same: uncertain ratings");
            case 4
                eval(['plot(SOA, dPar' num2str(i) '(9,:), ''b-*'')']); % old rating: same position
                title("Par " + i + " - diff side: joined ratings");
            case 5
                eval(['plot(SOA, dPar' num2str(i) '(7,:), ''b-*'')']); % old rating: same position
                title("Par " + i + " - diff: certain ratings");
            case 6
                eval(['plot(SOA, dPar' num2str(i) '(8,:), ''b-*'')']); % old rating: same position
                title("Par " + i + " - diff: uncertain ratings");
        end 
        xticks(SOA);
        xtickangle(90);
        xlabel('SOA');
        xlim([-350 350]);
        ylim([-0.05 1.1]);
        ylabel('mean ratings %');
    end
    
    fig_name = "confidence_plot_p" + i;
    png_name = fig_name + ".png";
    saveas(f, fig_name);
    saveas(f, png_name);
end

