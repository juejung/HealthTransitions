
// ----------------------------------------------------------------------------
//                      Path definitions
// ----------------------------------------------------------------------------
// Figure out the computer-hostname. Then set paths accordingly.
global jjOS "`c(hostname)'"
global jjuserName "`c(username)'"

disp "Run on Linux ..."
global MEPSdataDir   "/home/$jjuserName/Dropbox/AAA/Data/MEPS/MEPS_Annuals"
global HRSdataDir    "/home/$jjuserName/Dropbox/AAA/Data/HRS/Rand_1992_2016_V2"
global CPIdir        "/home/$jjuserName/Dropbox/AAA/Data/Indices"

global DatDir        "/home/$jjuserName/Dropbox/AAPapers/Jung/HealthMarkov/Stata/Data"
global StataOutDir   "/home/$jjuserName/Dropbox/AAPapers/Jung/HealthMarkov/Stata/Output"
global DoDir         "/home/$jjuserName/Dropbox/AAPapers/Jung/HealthMarkov/Stata/Do"
global PythonDir     "/home/$jjuserName/Dropbox/AAPapers/Jung/HealthMarkov/Python"
global AnacondaDir   "/home/$jjuserName/anaconda3/bin"
global TexDir        "/home/$jjuserName/Dropbox/AAPapers/Jung/HealthMarkov/Tex/Graphs"

