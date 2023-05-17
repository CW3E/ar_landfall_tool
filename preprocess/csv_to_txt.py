import pandas as pd
import numpy as np

textpts_fname = '../data/AK/AK_coastal_points.csv'
df = pd.read_csv(textpts_fname)
df = df.drop(['Unnamed: 0', 'name'], axis=1)
df = df[['lat-round', 'lon-round']]
df = df.sort_values(by=['lat-round', 'lon-round'])
idx = (df['lat-round'] < 60.) & (df['lon-round'] < 140.)
df = df.loc[idx]
df
np.savetxt(r'../data/AK/latlon_coast.txt', df.values, fmt='%.2f')