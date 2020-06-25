######################################################################
#
#   artemide for python3
#           A.Vladimirov (22.06.2020)
#
######################################################################

import artemide
import numpy

def initialize(fileName):
    """
    Initialization of artemide

    Parameters
    ----------
    fileName : string
        The path to the constants-file

    Returns
    -------
    None.

    """
    
    import os.path
    if os.path.exists(fileName):
    
        artemide.harpy.initialize(fileName)
        if artemide.harpy.started:
            pass
        else:
            print("Welcome to harpy -- the python interface for artemide")
    else:
        raise FileNotFoundError('consts-file '+fileName+'NOT FOUND')

def ShowStatistics():
    """
    Print the statistics

    Returns
    -------
    None.

    """
    artemide.harpy.showstatistics()

def setNPparameters(l):
    """
    Set NP parameters 

    Parameters
    ----------
    l : float array
        list of NP parameters

    Returns
    -------
    None.

    """
    if isinstance(l,list) or isinstance(l,numpy.ndarray):
        artemide.harpy.setlambda_main(numpy.asfortranarray(l))
    else:
        raise TypeError()
        
def setNPparameters_TMDR(l):
    """
    Setting NP parameters for the model of TMDR

    Parameters
    ----------
    l : float array or integer
        list of NP parameters, or the number of replica (if supported by model)

    Returns
    -------
    None.

    """
    if isinstance(l,list) or isinstance(l,numpy.ndarray):
        artemide.harpy.setlambda_tmdr(numpy.asfortranarray(l))
    elif isinstance(l, int):
        artemide.harpy.setreplica_tmdr(l)
    else:
        raise TypeError()
        

def setNPparameters_uTMDPDF(l):
    """
    Setting NP parameters for the model of uTMDPDF

    Parameters
    ----------
    l : float array or integer
        list of NP parameters, or the number of replica (if supported by model)

    Returns
    -------
    None.

    """
    if isinstance(l,list) or isinstance(l,numpy.ndarray):
        artemide.harpy.setlambda_utmdpdf(numpy.asfortranarray(l))
    elif isinstance(l, int):
        artemide.harpy.setreplica_utmdpdf(l)
    else:
        raise TypeError()
        
        
def setNPparameters_uTMDFF(l):
    """
    Setting NP parameters for the model of uTMDFF

    Parameters
    ----------
    l : float array or integer
        list of NP parameters, or the number of replica (if supported by model)

    Returns
    -------
    None.

    """
    if isinstance(l,list) or isinstance(l,numpy.ndarray):
        artemide.harpy.setlambda_utmdff(numpy.asfortranarray(l))
    elif isinstance(l, int):
        artemide.harpy.setreplica_utmdff(l)
    else:
        raise TypeError()
        
def setNPparameters_lpTMDPDF(l):
    """
    Setting NP parameters for the model of lpTMDPDF

    Parameters
    ----------
    l : float array or integer
        list of NP parameters, or the number of replica (if supported by model)

    Returns
    -------
    None.

    """
    if isinstance(l,list) or isinstance(l,numpy.ndarray):
        artemide.harpy.setlambda_lptmdpdf(numpy.asfortranarray(l))
    elif isinstance(l, int):
        artemide.harpy.setreplica_lptmdpdf(l)
    else:
        raise TypeError()

def setNPparameters_SiversTMDPDF(l):
    """
    Setting NP parameters for the model of SiversTMDPDF

    Parameters
    ----------
    l : float array or integer
        list of NP parameters, or the number of replica (if supported by model)

    Returns
    -------
    None.

    """
    if isinstance(l,list) or isinstance(l,numpy.ndarray):
        artemide.harpy.setlambda_siverstmdpdf(numpy.asfortranarray(l))
    elif isinstance(l, int):
        artemide.harpy.setreplica_siverstmdpdf(l)
    else:
        raise TypeError()


def varyScales(c1,c2,c3,c4):
    """
    Set new scale variation parameters

    Parameters
    ----------
    c1 : float
        Scale variation parameter c1
    c2 : float
        Scale variation parameter c2
    c3 : float
        Scale variation parameter c3
    c4 : float
        Scale variation parameter c4

    Returns
    -------
    None.

    """
    if not isinstance(c1,float):
        raise TypeError("c1 is not float")
    if not isinstance(c2,float):
        raise TypeError("c2 is not float")
    if not isinstance(c3,float):
        raise TypeError("c3 is not float")
    if not isinstance(c4,float):
        raise TypeError("c4 is not float")
        
    artemide.harpy.setscalevariation(c1,c2,c3,c4)

def _IsKinematicProper(s,qT,Q,y):
    """ Checks the point for the proper kinematics
    Especially for the correct domain of X.
    """
    gridX=0.00001
    if qT[0]>qT[1]:
        print('Wrong order of qT')
        return False
    if Q[0]>Q[1]:
        print('Wrong order of Q')
        return False
    if y[0]>y[1]:
        print('Wrong order of y')
        return False
    if qT[1]>Q[0]:
        print('qT (',qT[1],') > Q(', Q[0],')')
        return False
    if Q[1]>s:
        print('Q (',Q[1],') > s(', s,')')
        return False
    
    x1x2=(Q[0]**2+qT[0]**2)/s
    ymax=-numpy.log(numpy.sqrt(x1x2))
    ymin=-ymax
    
    if y[1]<ymin or y[0]>ymax:
        print('y ',y, 'is outside physical region ',[ymin,ymax])
        return False
    
    if y[1]>ymax:
        if x1x2<gridX:
            print('x outside of the grid')
            return False
    else:
        if numpy.sqrt(x1x2)*numpy.exp(-y[1])<gridX:
            print('x outside of the grid')
            return False
    
    if y[0]<ymin:
        if x1x2<gridX:
            print('x outside of the grid')
            return False
    else:
        if numpy.sqrt(x1x2)*numpy.exp(y[0])<gridX:
            print('x outside of the grid')
            return False
     
    x1x2=(Q[0]**2+qT[0]**2)/s
    ymax=-numpy.log(numpy.sqrt(x1x2))
    ymin=-ymax
    
    if y[1]<ymin or y[0]>ymax:
        print('y ',y, 'is outside physical region ',[ymin,ymax])
        return False
    
    if y[1]>ymax:
        if x1x2<gridX:
            print('x outside of the grid')
            return False
    else:
        if numpy.sqrt(x1x2)*numpy.exp(-y[1])<gridX:
            print('x outside of the grid')
            return False
    
    if y[0]<ymin:
        if x1x2<gridX:
            print('x outside of the grid')
            return False
    else:
        if numpy.sqrt(x1x2)*numpy.exp(y[0])<gridX:
            print('x outside of the grid')
            return False
    
    return True

def setPDFreplica(n):
    """
    Changes the replica for PDF input.
    
    This is a temporary function will be changed in future versions

    Parameters
    ----------
    n : Integer
        Number of PDF replica

    Returns
    -------
    None.

    """
    if not isinstance(n,int):
        raise TypeError()
    artemide.harpy.setpdfreplica(n)
   
###############################################################################
class DY:
        """Static class for evaluation of DY cross-section
        """

        @staticmethod
        def xSec(process,s,qT,Q,y,includeCuts,CutParameters=None,Num=None):
                """Cross-section for DY integrated over bin
                      
                Arguments: (process,s,qT,Q,y,includeCuts,CutParameters=None,Num=4)
                process         = (int, int, int) (see definition in artemide manual)
                s               = Mandelshtam variable s
                qT              = (qT-Min,qT-Max) boundaries of bin in qT
                Q               = (Q-Min,Q-Max) boundaries of bin in Q
                y               = (y-Min,y-Max) boundaries of bin in y
                includeCuts     = True/False to include leptonic cuts
                CutParameters   = (real,real,real,real) must be if includeCuts=True (see definition in artemide manual)
                Num             = even integer, number of section of qt-integration (defaul=4)
                """

                
                if not includeCuts:
                        cc=[0,0,0,0]
                else:
                        cc=CutParameters
                
                if Num==None:
                    return artemide.harpy.dy_xsec_single(\
                        numpy.asfortranarray(process),\
                        s,\
                        numpy.asfortranarray(qT),\
                        numpy.asfortranarray(Q),\
                        numpy.asfortranarray(y),\
                        includeCuts,\
                        numpy.asfortranarray(cc))
                else:
                    return artemide.harpy.dy_xsec_single(\
                        numpy.asfortranarray(process),\
                        s,\
                        numpy.asfortranarray(qT),\
                        numpy.asfortranarray(Q),\
                        numpy.asfortranarray(y),\
                        includeCuts,\
                        numpy.asfortranarray(cc),\
                        Num)
                    
            
        @staticmethod
        def xSecList(process,s,qT,Q,y,includeCuts,CutParameters):
            return artemide.harpy.dy_xsec_list(numpy.asfortranarray(process),\
                                               numpy.asfortranarray(s),\
                                               numpy.asfortranarray(qT),\
                                               numpy.asfortranarray(Q),\
                                               numpy.asfortranarray(y),\
                                               numpy.asfortranarray(includeCuts),\
                                               numpy.asfortranarray(CutParameters),\
                                               len(s))
        @staticmethod
        def xSecListBINLESS(process,s,qT,Q,y,includeCuts,CutParameters):
            """ The evaluation of cross-section at a single point. 
            """
            #print len(s)
            
            return artemide.harpy.dy_xsec_binless_list(numpy.asfortranarray(process),\
                                               numpy.asfortranarray(s),\
                                               numpy.asfortranarray(qT),\
                                               numpy.asfortranarray(Q),\
                                               numpy.asfortranarray(y),\
                                               numpy.asfortranarray(includeCuts),\
                                               numpy.asfortranarray(CutParameters),\
                                               len(s))
                                
###############################################################################
class SIDIS:
        """Static class for evaluation of SIDIS cross-section
        """

        @staticmethod
        def xSec(process,s,pT,z,x,Q,includeCuts,CutParameters=None,masses=[0.938,0.130]):
                """Cross-section for DY integrated over bin
                      
                Arguments: (process,s,qT,z,x,Q,includeCuts,CutParameters=None)
                process         = (int, int, int) (see definition in artemide manual)
                s               = Mandelshtan variable s
                pT              = (pT-Min,pT-Max) boundaries of bin in pT
                z               = (z-Min,z-Max) boundaries of bin in z
                x               = (x-Min,x-Max) boundaries of bin in x
                Q               = (Q-Min,Q-Max) boundaries of bin in Q                
                includeCuts     = True/False to include leptonic cuts
                CutParameters   = (real,real,real,real) must be if includeCuts=True (see definition in artemide manual)
                """
                
                if not includeCuts:
                        cc=[0,0,0,100]
                else:
                        cc=CutParameters
                
                
                return artemide.harpy.sidis_xsec_single_withmasses(\
                    numpy.asfortranarray(process),\
                    s,\
                    numpy.asfortranarray(pT),\
                    numpy.asfortranarray(z),\
                    numpy.asfortranarray(x),\
                    numpy.asfortranarray(Q),\
                    includeCuts,\
                    numpy.asfortranarray(cc),
                    numpy.asfortranarray(masses))
                    
            
        @staticmethod
        def xSecList(process,s,pT,z,x,Q,includeCuts,CutParameters,masses):  
            return artemide.harpy.sidis_xsec_list(numpy.asfortranarray(process),\
                                               numpy.asfortranarray(s),\
                                               numpy.asfortranarray(pT),\
                                               numpy.asfortranarray(z),\
                                               numpy.asfortranarray(x),\
                                               numpy.asfortranarray(Q),\
                                               numpy.asfortranarray(includeCuts),\
                                               numpy.asfortranarray(CutParameters),\
                                               numpy.asfortranarray(masses),\
                                               len(s))
            
        @staticmethod
        def xSecListBINLESS(process,s,pT,z,x,Q,masses):  
            """ The evaluation of cross-section at a single point. 
            
            No binning effects.
            Consiquetly, no cuts.
            """
            return artemide.harpy.sidis_xsec_binless_list(numpy.asfortranarray(process),\
                                               numpy.asfortranarray(s),\
                                               numpy.asfortranarray(pT),\
                                               numpy.asfortranarray(z),\
                                               numpy.asfortranarray(x),\
                                               numpy.asfortranarray(Q),\
                                               numpy.asfortranarray(masses),\
                                               len(s))
                                
