import pandas as pd

params_file='/home/cnari/cnari_projects/sports_2014/Scripts/preproc/sports_preproc_params.csv'
params_frame=pd.read_csv(params_file, index_col='sub')
params_frame.save('/home/cnari/cnari_projects/sports_2014/Scripts/preproc/sports_preproc_params.pck')