function create_all()

%sequences = {'00','01','02','03','05','06','07','08','09','10'};
sequences = {'05','06','07','08','09','10'};

DATASET_DIR = '/home/kreimer/KITTI/';

parfor i=1:length(sequences)
    feature.harris_tri_match(DATASET_DIR, sequences)
end