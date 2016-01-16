function updated_tracks = read_tracks(filename, tracks)

fd = fopen(filename, 'r');
A = fscanf(fd,'%d %d %d %d %d', [5 Inf])';
fclose(fd);

updated_tracks = cell(0);
for i = 1:size(A,1)
    cur_ind = A(i,1);
    
    found = false;
    for j=1:size(tracks,2)
        ind = tracks{2,j}==cur_ind;
        
        if any(ind)
            if size(updated_tracks,2) < j+1
                updated_tracks{1,j+1} = [];
            end
            updated_tracks{1,j+1}(:,end+1) = [tracks{1,j}(:,ind); A(i,2:end)'];
            updated_tracks{2,j+1}(1,end+1) = tracks{2,j}(1,ind);
            found = true;
            break;
        end
    end
    if ~found
        if isempty(updated_tracks)
            updated_tracks{1,1}(:,1) = A(i,2:end)';
            updated_tracks{2,1}(:,1) = A(i,1);
        else
            updated_tracks{1,1}(:,end+1) = A(i,2:end)';
            updated_tracks{2,1}(:,end+1) = A(i,1);            
        end
    end
end
end