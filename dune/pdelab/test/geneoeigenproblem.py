import scipy
import scipy.io
import scipy.linalg
import numpy
import matplotlib.pyplot as plt
import cmath

FineMatrixExterior = scipy.sparse.coo_matrix.todense(scipy.io.mmread("0FineMatrixExterior.mm"))
PartUnityOverlapMatrixPartUnity = scipy.sparse.coo_matrix.todense(scipy.io.mmread("0PartUnityOverlapMatrixPartUnity.mm"))
PartUnityExteriorMatrixPartUnity = scipy.sparse.coo_matrix.todense(scipy.io.mmread("0PartUnityExteriorMatrixPartUnity.mm"))
DTNFineMatrix = scipy.sparse.coo_matrix.todense(scipy.io.mmread("0DTNFineMatrix.mm"))
#FineMatrixExterior = scipy.io.mmread("0FineMatrixExterior.mm")
#PartUnityOverlapMatrixPartUnity = scipy.io.mmread("0PartUnityOverlapMatrixPartUnity.mm")

#print FineMatrixExterior
#[w, vl] = scipy.linalg.eig(FineMatrixExterior, PartUnityOverlapMatrixPartUnity)
#print w
#print vl


plt.ylabel('Eigenvalues')
plt.yscale('log')
plt.title('Spectra of eigenproblems')

def calcEigenproblem(A, B, title):
  print("Shape:")
  print(A.shape)
  #[w, vl] = scipy.linalg.eig(FineMatrixExterior)
  [w, vl] = scipy.linalg.eig(a=A, b=B)
  #plt.figure()
  #plt.plot(numpy.sort(w)[::-1], label = title)
  plt.plot(numpy.sort(w.real), label = title)

  print (w.shape)
  print(numpy.sort(w.real))

calcEigenproblem(FineMatrixExterior, PartUnityOverlapMatrixPartUnity, "GenEO")
calcEigenproblem(FineMatrixExterior, PartUnityExteriorMatrixPartUnity, "GenEO with interior")
calcEigenproblem(FineMatrixExterior, DTNFineMatrix, "Dirichlet to Neumann")

plt.legend()
plt.draw()
plt.show()
#scipy.sparse.linalg.eigs(A = FineMatrixExterior, M = PartUnityOverlapMatrixPartUnity)
