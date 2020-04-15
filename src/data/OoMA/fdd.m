function [SV,F,Phi,I,SV_nums] = fdd(Y,win,fs,peak_freqs)
% [SV,F,Phi,I,SV_nums] = FDD(Y,win,fs,peak_freqs)
%   Frequency Domain Decomposition (FDD) user interface for peak-picking.
%
%   INPUTS:
%   Y           time-domain sensor data array
%   win         (optional) window to use for computing the CMIF, use [] if
%               unspecified
%   fs          sampling frequency used for the data Y
%   peak_freqs  (optional) frequencies of the peaks in the power spectrums
%               or CMIF plots to use for peak peaking rather than
%               interactive mode, use [] if unspecified (manual mode)
%
%   OUTPUTS:
%   SV          singular-value components from the CMIF
%   F           frequency vector corresponding to the singular values SV
%   Phi         matrix of mode shape vectors, where Phi(1,:) is the first
%               mode shape vector
%   I           indices of the frequency vector F corresponding to the
%               identified modes in Phi
%   SV_nums     the CMIF component number (singular value vector)
%               corresponding to the mode shape in Phi, where Phi(i,:)
%               corresponds to SV_nums(i)
%
% REFERENCES:
% [1]   Brincker, Rune, Lingmi Zhang, and P. Andersen. "Modal
%       identification from ambient responses using frequency domain
%       decomposition." Proc. of the 18th International Modal Analysis
%       Conference (IMAC), San Antonio, Texas. 2000.

% generate the CMIF
[SV,F,U] = cmif(Y,win,1/fs);
h = figure();

% only show up to 4 singular value vectors from the CMIF
if size(SV,1) >= 4
    n_components = 4;
else
    n_components = size(SV,1);
end

% default select the first singular value vector
sel_component = 1;

% set up the UI
axes
pos = get(gca,'position');
set(gca,'position',[pos(1) pos(2) .6 pos(3)])
uicontrol('units','normalized','style','pushbutton','string','Select Peak','position',[.75 .65 .2 .1],'callback',@setpeak)
uicontrol('units','normalized','style','pushbutton','string','Clear Peak','position',[.75 .55 .2 .1],'callback',@clearpeak)
uicontrol('units','normalized','style','pushbutton','string','Finish','position',[.75 .2 .2 .1],'callback',@finish)
uicontrol('units','normalized','style','text','string','CMIF Component','position',[.75 .85 .2 .05])
uicontrol('units','normalized','style','popup','string',{1:n_components},'position',[.75 .75 .2 .1],'callback',@component)
uicontrol('units','normalized','style','text','string','Status: ','position',[.75 .45 .2 .05],'horizontalalignment','left')
status_text = uicontrol('units','normalized','style','text','string','','position',[.75 .35 .2 .1],'horizontalalignment','left');

% initialized peak picking arrays
peak_vals = [];
I = [];
SV_nums = [];

% draw the plot area
refreshPlot

% wait until the user presses the finish button or closes the figure
waitfor(h)

% create the mode shape matrix
Phi = zeros(size(SV,1),length(I));
[~,sort_idx] = sort(I);
I = I(sort_idx);
SV_nums = SV_nums(sort_idx);
for idx = 1:length(I)
    Phi(:,idx) = U(:,SV_nums(idx),I(idx));
end

    function setpeak(source,event)
        status_text.String = 'Drag rectangle to select a peak';
        rect_pos = getrect(h);
        x1 = rect_pos(1);
        x2 = x1 + rect_pos(3);
        idx1 = find(x1>F,1,'last');
        idx2 = find(F>x2,1,'first');
        [pk,pk_idx] = max(SV(sel_component,idx1:idx2));
        peak_vals = [peak_vals pk];
        pk_idx = pk_idx + idx1 - 1;
        I = [I pk_idx];
        SV_nums = [SV_nums sel_component];
        status_text.String = ['Peak selected at ' num2str(round(F(pk_idx),2)) ' Hz'];
        refreshPlot
    end

    function clearpeak(source,event)
        if ~isempty(I)
            I(end) = [];
            SV_nums(end) = [];
            peak_vals(end) = [];
            status_text.String = 'Last peak deleted';
        end
        refreshPlot
    end

    function finish(source,event)
        % just close the window
        close(h)
    end

    function component(source,event)
        sel_component = source.Value;
        status_text.String = ['CMIF Component ' num2str(sel_component) ' selected'];
        refreshPlot
    end

    function refreshPlot
        fdd_lines = plot(F,10*log10(abs(SV(1:n_components,:))));
        line_colors = lines(2);
        hold on
        for i = 1:n_components
            if i == sel_component
                fdd_lines(i).Color = line_colors(1,:);
            else
                fdd_lines(i).Color = line_colors(2,:);
            end
        end
        plot(F(I),10*log10(peak_vals),'ko');
        hold off
        xlim([F(1) F(end)])
        xlabel('Frequency (Hz)')
        ylabel('CMIF Magnitude (dB)')
        title('Frequency Domain Decomposition')
    end


end