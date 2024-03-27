import rpy2.robjects as robjects
import pandas as pd

# Load the RDS file
rds_file_path = "/Users/amaros/Documents/mgCSTs2/RDS_files/vog.list.mgss.pa.RDS"

# Load the RDS file into an R object
r_data = robjects.r['readRDS'](rds_file_path)

# Convert the R object to a Python object
py_data = robjects.conversion.rpy2py(r_data)

# Now you can work with `py_data` in Python
df = pd.DataFrame(py_data[1])
print(df.keys())
#df.to_csv("test.csv", index=False)

