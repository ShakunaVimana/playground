class FirstClass:
  def setdata(self,value):
    self.data = value
  def display(self):
    print(self.data)

x = FirstClass()
x.setdata(3.14)
x.display()
x.data = "King Arthur"
x.display()

class SecondClass(FirstClass):
  def display(self): # Changes display
    print('Current value = "%s"' % self.data)

class ThirdClass(SecondClass): # Inherit from SecondClass
  def __init__(self, value): # On "ThirdClass(value)"
    self.data = value
  def __add__(self, other): # On "self + other" #change define of "+"
    return ThirdClass(self.data + other)
  def __str__(self): # On "print(self)", "str()" #change define of "str"
    return '[ThirdClass: %s]' % self.data 
  def mul(self, other): # In-place change: named
    self.data *= other

a = ThirdClass('abc') # __init__ called
a.display() # Inherited method called
print(a) # __str__: returns display string
b = a + 'xyz' # __add__: makes a new instance
b.display() # b has all ThirdClass methods
print(b) # __str__: returns display string
a.mul(3) # mul: changes instance in place
print(a)
