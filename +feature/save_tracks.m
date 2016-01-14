function save_tracks(filename, tracks, x1, x2, m12, mt)

fd = fopen(filename, 'w+');
for i=1:size(tracks,2)
    for j=1:size(tracks{1,i},2)
        f1 = tracks{1,i}(end,j);
        if i>1
            ind = find(mt(1,:) == f1,1);
            f2 = mt(4,ind);
        else
            ind = find(m12(1,:) == f1,1);
            f2  = m12(2,ind);
        end
        
        c1 = x1(:,f1);
        c2 = x2(:,f2);
        fprintf(fd,'%d %d %d %d %d\n', tracks{2,i}(j), c1(1), c1(2), c2(1), c2(2));
    end
end
fclose(fd);
end
