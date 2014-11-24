import numpy as np
import pylab as pl
from math import exp
class NumericalMethods(object):
	"""
	A class comparing three numerical methods:
	forward Euler, backward Euler and Cranck
	Nickolsson scheme. The theta scheme is 
	implemented: equation: u'(t) = -au(t) for a>0
				u(0) = I
	"""

	def __init__(self,initialValue,timeStep,a,domainSize, theta = None):

		self.initialValue = initialValue   # I
		self.timeStep         = timeStep #dt
		self.a                = a # constant a in ode
		self.domainSize       = domainSize # T
		self.theta            = theta      # theta in the tetha rule

		self. N = int(self.domainSize/self.timeStep)
		self. norm = 0
		self.error = 0
		self.u = np.zeros(self.N+1)
		self.ue = np.zeros(self.N+1)
		self.u[0] = self.initialValue
	
	def forwardEuler(self):
		for n in range(0,self.N):
			self.u[n+1] = self.u[n]*(1-self.a*self.timeStep)
			
			
	def backwardEuler(self):
		for n in range(0,self.N):
			self.u[n+1] = (self.u[n])/(1+self.a*self.timeStep)		

	def crankNikolson(self):
		for n in range(0,self.N):
			self.u[n+1] = self.u[n]*( (1-0.5*self.a*self.timeStep)/(1+0.5*self.a*self.timeStep)   )
			
	def exactSolution(self):
		t = np.linspace(0,self.domainSize,self.N)
		self.ue[0:self.N] = self.initialValue*np.exp(-self.a*t[0:self.N])

 	def thetaRule(self):
 		for n in range(0,self.N):
 			self.u[n+1] = (   (1-(1-self.theta)*self.a*self.timeStep)/(1+self.theta*self.a*self.timeStep) )*self.u[n]


	def verifications(self,s):
		error = np.linalg.norm(self.ue-self.u, ord=2)
		print s,error

	def plotSolution(self,s):
		"""
		customize plot for each method
		"""
 		t = np.linspace(0,self.domainSize,self.N+1)
 		pl.plot(t,self.u,t,self.ue)
 		pl.title(s)
 		pl.ylim([-2,self.initialValue+0.5])
 		pl.show()
		


def main():
	"""
	initialization
	"""
	I, dt, a,T = 3.5,0.25,1,10
	
 	numericalmethod= NumericalMethods(I,dt,a,T)
 	
 	
 	numericalmethod.exactSolution() # call method implemeting exactsolution
 	
 	numericalmethod.forwardEuler()
 	numericalmethod.plotSolution("ForwardEuler ,dt = %f"%dt)
 	
 	
 	numericalmethod.backwardEuler()
 	numericalmethod.plotSolution("backward Euler, dt = %f"%dt)
 	
 	numericalmethod.crankNikolson()
 	numericalmethod.verifications("L2 norm for crankNikolson: ")
 	numericalmethod.plotSolution("crankNikolson, dt = %f"%dt)
 	
 	ob = NumericalMethods(I,dt,a,T,0)
 	ob.thetaRule()
 	ob.plotSolution(" theta rule for theta = 0")
 	
 	
 	#crankNikolson.plotSolution()
if __name__=="__main__":
	main()
