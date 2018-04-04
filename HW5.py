import numpy
import math

def Rom(f,a,b,n):
  r = numpy.array([[0]*(n+1)] * (n+1),float)
  h = b - a
  r[0,0] = 0.5 * h * (f(a)+f(b))

  powerOf2 = 1
  for i in range(1, n+1 ):
      h = 0.5 * h
      sum = 0.0
      powerOf2 = 2 * powerOf2
      for k in range(1,powerOf2,2):         
        sum = sum + f(a+k*h)
        r[i,0] = 0.5 * r[i-1,0] + sum * h
        powerOf4 = 1
        for j in range(1,i+1):
            powerOf4 = 4 * powerOf4
            r[i,j] = r[i,j-1] + (r[i,j-1]-r[i-1,j-1])/(powerOf4 - 1)
  print("The approximate value of the integral is", r[n-1][n-1])


  
def Simp(f,a,b,n):
    h = float(b-a)/n
    odd = 0.0
    even = 0.0
    for i in range(1, n, 2):
        xk = a + i*h
        odd += f(xk)
    for i in range(2, n, 2):
        xk = a + i*h
        even += f(xk)
    answer = 2*even + 4*odd + f(a) + f(b)
    print("The approximate value of the integral with", n, "subintervals is", (h/3)*answer)

def Gauss3(f,a,b):
  Matrix = [[-math.sqrt(3/5),5/9],[0, 8/9], [math.sqrt(3/5),5/9]]
  one = Matrix[0][1] * f(((b-a)*Matrix[0][0] + (a+b))/2)
  two = Matrix[1][1] * f((a+b)/2)
  three = Matrix[2][1] * f(((b-a)*Matrix[2][0] + (a+b))/2)
  y = (one + two + three) * (b-a)/2
  print(y)

    
def function(x): return math.tan(x**3)

def numfour(x): return math.exp(-x**2)

def numfive(x): return (10**-4)/((x-(math.pi/2))**2 + (10**-8))

def numsix(x): return (x**2 - 1)**(1./3.)/ math.sqrt(math.sin(math.exp(x)-1))

print("Number Three Test:")
Rom(function,0,1,4)
Gauss3(function, 0,1) 
print("")

print("Number Four:")
Simp(numfour, 0,1,16)
Rom(numfour, 0,1,4)
#With the simpsons rule the error comes out to -1.244 * 10^-7,and the Romberg comes out to 1.14510 * 10^-7. Both give really good approximations, however the Romberg gives a slightly better answer. 

print("")
print("Number Five:")
Simp(numfour, 0,2,16)
Simp(numfive, 0,2,128)
Simp(numfive, 0,2,1024)
Rom(numfive, 0,2,4)
Rom(numfive, 0,2,7)
Rom(numfive, 0,2,10)
print("")
Simp(numfive, 0,2,32768)
Rom(numfive, 0,2,15)
#The more subintervals you use with Simpsons and Romberg, the more accurate the approximation becomes. Mainly due to the spike in the graph between 1 and 2. 

print("")
print("Number Six:")
Gauss3(numsix,0,1)
#Division by zero becomes a problem when trying to use the trapezoid or simpsons rule. 
