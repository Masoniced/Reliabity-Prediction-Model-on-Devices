import pandas as pd

file_name = r"C:\Users\Sen\Desktop\Raw-Curves Files\Logan-7C_w12_TDDB_25C_Compiled Raw.txt"
File = pd.read_csv(file_name, sep ='\t', header = 0)
Columns = File.columns
Result = pd.DataFrame(columns = Columns)
Ini_key = File.iat[0,2]
criteria = 0
Low_resistance = 1000

for i in range(1, len(File.index)):
	if File.iat[i,2] != Ini_key:
		Ini_key = File.iat[i,2]
		criteria = 0
	else:
		if criteria == 0 and File.iat[i,11] < Low_resistance:
			Result.loc[Result.shape[0]] = File.iloc[i]
			criteria = 1
		else:
			continue

pd.DataFrame(Result).to_csv('Low_resistance'+'.txt', index=False)

