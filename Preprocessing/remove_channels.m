function [chans_to_remove] = remove_channels(pk,labels)
        chans_to_remove = [];
        for chani = 1:length(pk.chanlocs)
            if isempty(cell2mat(strfind(labels,pk.chanlocs(chani).labels)))
                chans_to_remove = [chans_to_remove,chani];        
            end
        end
    end