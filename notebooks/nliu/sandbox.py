# %%
import json
import numpy as np

# %%
with open("./../../mira/dkg/resources/probonto.json", "r") as f:
    probonto = json.load(f)

# %%
ciemss = []
with open("./test", "r") as f:
    for l in f:
        ciemss.append(l)

ciemss = np.array(ciemss)

ciemss = ciemss[0::3]
print(len(ciemss))

# %%

