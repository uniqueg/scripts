#!/usr/bin/env python

# Alexander Kanitz
# 04-AUG-2016

import sys

n = int(sys.argv[1])

def prime_factors(n):
  fact = []
  rem = 0 
  for i in range(2,n+1):
    if not n % i:
      fact.append(i)
      rem = n // i
      if rem > 1:
        fact = fact + prime_factors(rem)
        return(fact)
  return(fact)

def add_missing_factors(x,y):
  a = x[:]
  b = y[:]
  union = []
  for i in list(b):
    if i in a:
      a.remove(i)
    else:
      union.append(i)
    b.remove(i)
  union = x + union
  union.sort()
  return union

def div_by_all_int_up_to_n(n):
  minfact = []
  fact = []
  for i in range(2, n+1):
    fact = prime_factors(i)
    minfact = add_missing_factors(minfact, fact)
  return(minfact)

def main():
  fact = div_by_all_int_up_to_n(n)
  prod = reduce(lambda x, y: x*y, fact)
  print 'Number:'
  print n
  print 'Factors:'
  print fact
  print 'Product:'
  print prod

main()
