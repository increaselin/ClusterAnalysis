__Readme.first.txt

 These scripts and programs are used to correct equivalent width measures of the oxygen triplet for non-local thermal equilibrium effects.

 The scripts will run in a "batch" mode on all paired .log and .m files in the directory.
 
 The original paper (Takeda 2003) can be viewed at:
 http://www.aanda.org/index.php?option=com_article&access=standard&Itemid=129&url=/articles/aa/full/2003/16/aa2727/aa2727_online.html
 
 Usage:
 ./addNLTEOEffects.py

Input files:
<starname>.log : A file containing any number of line measurements, in the output form of the "splot" function in IRAF.
<starname>.m : A corresponding atmospheric model file from the "mspawn" function. The temperature, logg and microturbulent velocity parameters are used by these scripts.

Output file:
<starname>_NLTEO.log : A modified log file, containing the same line width measures as the original .log file, with only the oxygen triplet line widths changed. Note: The scripts will also add a notation to the first comment line of the log file to indicate that the changes have been made.

Data files:
__Readme.first.txt : This file
NLTEOxyCorrectionDB : A MySQL database containing the correction tables from (Takeda 2003).
Takeda_Oxy_Table3.tab : The text data from the (Takeda 2003) paper. This file is not actually used in the processing, but is included in case the user wishes to make adjustments in the future.
createODB.py : A python script used to construct the database of corrections (NLTEOxyCorrectionDB) from a text table. Use if there are changes to the NLTE adjustments from those in (Takeda 2003).