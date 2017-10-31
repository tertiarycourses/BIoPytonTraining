file=open('C:\\Users\\user\\Desktop\\Python\\pytTUTORIALS\\datasets\\hitler.txt', mode='r', encoding="utf8")


file.read()


file.closed   # check if our file is closed > returns false



file.close()

file.closed  # now returns true









########## reading csv files

import pandas as pd



heart='C://Users//user//Desktop//tertiary//PYTprogramming//BioPython//datasets//titanic.csv'


df=pd.read_csv(heart)

df.head()





### slicing

heart2=pd.read_csv(heart, nrows=5, header=None)

  # reading only 5 rows

myheart=heart[   :,3:4]   # selecting columns 3 to 4

myheart2=heart[['age','trestbps','chol','thalach']]


