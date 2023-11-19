import os

import pandas as pd


def find_fonts():
    for filename in os.listdir('fonts'):
        df = pd.read_csv(f"fonts/{filename}")
        if len(df) > 5_000:
            print(filename.strip('.csv'), len(df))


def combine_fonts():
    result_data = []

    fonts = ['ARIAL', 'TIMES', 'CALIBRI', 'MONEY', 'NUMERICS']
    for font in fonts:
        df = pd.read_csv(f"fonts/{font}.csv")
        df_selected = df.iloc[:, 12:412]
        df_selected['font'] = font
        result_data.append(df_selected.sample(frac=1, random_state=42)[:10_000])
    final_df = pd.concat(result_data, ignore_index=True)
    shuffled = final_df.sample(frac=1, random_state=42)
    return shuffled[['font'] + [col for col in shuffled.columns if col != 'font']]


# find_fonts()
df = combine_fonts()
df.to_csv('fonts.csv', index=False)
