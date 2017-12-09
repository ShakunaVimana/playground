class Person:
 def __init__(self, name, job):
   self.name = name
   self.job = job
 def info(self):
   return (self.name, self.job)

rec1 = Person('mel', 'trainer')
rec2 = Person('vls', 'developer')

print rec1.job,rec2.info()

