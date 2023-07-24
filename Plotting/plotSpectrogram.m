function plotSpectrogram(c,t,f,varargin)

params = parseinputs(varargin{:});

surf(t,f,c,EdgeColor='none');
view(0,90);
shading interp;colormap(parula(128));
axis tight;
h = colorbar;
h.Label.String = 'Power in dB';
if isempty(params.xlab) && isempty(params.ylab)
    ylabel('Frequency (Hz)');xlabel('Time (s)');
else
    xlabel(params.xlab); ylabel(params.ylab);
end
if ~isempty(params.PlotTitle)
    title(params.PlotTitle);
else
    title('Wavelet based Spectogram');
end

%----------------------------------------------------------------
function params = parseinputs(varargin)
    
    params.PlotTitle = [];
    params.xlab = [];
    params.ylab = [];

    if isempty(varargin)
        return;
    end
    Len = length(varargin);
    if (Len==1)
        params.PlotTitle = varargin{1};
    end
    if (Len == 3)
        params.PlotTitle = varargin{1};
        params.xlab = varargin{2};
        params.ylab = varargin{3};
    end
       