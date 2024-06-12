# %%
import pandas as pd
import numpy as np
import os
import json
import matplotlib.pyplot as plt

# %%
DATA_PATH = r"https://raw.githubusercontent.com/ciemss/program-milestones/epi-scenario-1/18-month-milestone/hackathon/epi/Scenario%201%20Supplemental/"
dataset1 = os.path.join(DATA_PATH, "UK_compartments.csv")

# Load datasets 
d1 = pd.read_csv(dataset1, index_col=0)

# %%
grouped_sum = d1.groupby(['t', 'compartment', 'age_group'])['y'].sum()

reshaped_d1 = grouped_sum.unstack(fill_value=0)
reshaped_d2 = reshaped_d1.unstack(fill_value=0)

# %%
new_column_names = ['{}_{}'.format(compartment, age) if compartment else 't' for age, compartment in reshaped_d2.columns]
reshaped_d2.columns = new_column_names

reshaped_d2 = reshaped_d2.reset_index()


# # %%
# # Group by 't' and 'compartment', then sum 'y' values
# grouped = df.groupby(['t', 'compartment'])['y'].sum()


# # Reshape the DataFrame and fill missing values with 0
# df_new = grouped.unstack(fill_value=0)

# # Make a column 'I' which sums 'ICase' and 'IMild' and delete the original columns
# df_new['I'] = df_new['ICase'] + df_new['IMild']
# df_new.drop(columns=['ICase', 'IMild'], inplace=True)

# # Output to CSV
# df_new.to_csv('UK_SEIRD.csv')

# # Redefine the dataset
# dataset1 = 'UK_SEIRD.csv'