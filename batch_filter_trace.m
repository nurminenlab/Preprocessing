function filtr_data = batch_filter_trace(filename,ops,type,order,varargin)
% HP filters the data in 'filename' and writes the data to a new
% file
% based on Kilosort - preprocessData.m


if nargin == 6
    coff = varargin{1};
    l_coff = varargin{2};
else
    coff = varargin{1};
end

if exist([filename,'hp']) == 2

    fid        = fopen([filename,'hp'], 'r');
    filtr_data = fread(fid,[ops.NchanTOT Inf],'*int16');
    fclose(fid);

else
    
    NTbuff = ops.NT + 4*ops.ntbuff;
    ibatch = 0;
    fid         = fopen(filename, 'r');
    hpflter_fid = fopen([filename,'hp'], 'w');

    if type == 'high'
        [b1, a1] = butter(order, coff/ops.fs*2, type);
    elseif type == 'bandpass'
        [b1, a1] = butter(order, [l_coff/ops.fs*2, coff/ops.fs*2], type);
    elseif type == 'low'
        [b1, a1] = butter(order, coff/ops.fs*2, type);
    end
    
    while 1
        ibatch = ibatch + 1;
    
        offset = max(0, 2*ops.NchanTOT*NTbuff*(ibatch-1));
        fseek(fid, offset, 'bof');
        buff = fread(fid, [ops.NchanTOT NTbuff], '*int16');
        if isempty(buff)
            break;
        end

        dataRAW = gpuArray(buff);
        dataRAW = dataRAW';
        dataRAW = single(dataRAW);
        dataRAW = dataRAW(:, ops.chanMap);
        
        datr = filter(b1, a1, dataRAW);
        datr = flipud(datr);

        datr = filter(b1, a1, datr);
        datr = flipud(datr);
        
        fwrite(hpflter_fid, gather(datr)', 'int16');
    end
    
    fclose(fid);
    fclose(hpflter_fid);
    fid        = fopen([filename,'hp'], 'r');
    filtr_data = fread(fid,[ops.NchanTOT Inf],'*int16');

end
