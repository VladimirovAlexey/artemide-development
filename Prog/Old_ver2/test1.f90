!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	Elementary program that read old-constants-file (possibly for older version of artemide)
!!		and create a new (for current version of artemide)
!!
!!						A.Vladimirov (08.07.2019)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_setup
implicit none

CHARACTER(*),parameter::path='/home/vla18041/LinkData2/arTeMiDe_Repository/Constants-files/'
CHARACTER(*),parameter::fileIN='const-DYfit18_LO'
CHARACTER(*),parameter::fileOUT='const-DYfit18_LO+'

 call artemide_Setup_fromFile(fileIN,prefix=path)
 call CreateConstantsFile(fileOUT,prefix=path)
 
end program example