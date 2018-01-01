import scipy
import scipy.misc as sc
from math import factorial

sum = 0
for i in range(1,12):
  print scipy.misc.comb(12,i)*(i)**44
  sum = sum + factorial(12)*scipy.misc.comb(12,i)*(i)**44

x=  (12**56-sum)*1.0
y = 12**56
x/y
print x,y,x/y

z = factorial(56)*sc.comb(55,11)
print z*1.0/y

stirling = 0
for j in range(1,13):
	stirling = stirling + (-1)**(12-j)*sc.comb(12,j)*j**56

print stirling/y

