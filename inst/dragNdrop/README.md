This file is for admins/support-staff who want to configure and install the PTXQC package
and provide an automated drag'n'drop solution to their MaxQuant users (incl. themselves) for creating QC reports.


We recommend using a shared network drive as installation folder to which all users have read access. Alternatively you can copy the folder structure we
are about to a folder on your local machine - however, if you want to use the QC reporting from multiple PC's just use a network folder.

When you are done, provide the users with the [user_manual.pdf] from this folder (or write your own).
This will show them how to invoke the QC (it's very easy - see 'Usage' below).

### Installation
 - create clean R installation by downloading R from [http://www.r-project.org] and installing it. Make sure to install the x64 version for 64bit support.
 - start the R GUI and execute the following commands to install the `PTXQC` package (dependencies should be installed automatically):

        install.packages("devtools")
        library("devtools")             ## this might give a warning like 'WARNING: Rtools is required ...'. Ignore it.
        source("http://bioconductor.org/biocLite.R")
        biocLite("Biobase)
        install_github("cbielow/PTXQC") 
        help(package="PTXQC")           ## all done; check out the documentation


 - if R installed the new packages to a custom location (you can tell by the console output during package installation which would say something like: 'Installing package into ‘K:/R/win-library/3.1’')
   , move them to the .\library subfolder within your R installation folder (e.g., C:\Program Files\R\R-3.1.0\library - your path might differ a little)
 - copy the folder `QC-dragdrop` (from <R-installation-dir>\library\PTXQC\inst\dragNdrop) to your target location, usually some network or local drive (e.g. <target-drive> is Z:\proteomics\QC-dragdrop)
 - copy the whole R installation directory (e.g. c:\program files\R\R-3.1.0) into the `QC-dragdrop\\_internal` sub-folder
 - rename your newly copied R installation directory to R-3.1.0 (or edit the createQC_dragNdrop.bat to match your directory name)
 
Now, you should have the following structure

    <target-drive>\QC-dragdrop\
                               \createQC_dragNdrop.bat
                               \createQC_dragNdrop_withYAML.bat
                               \_internal\
                                         \R-3.1.0  (version may differ)
                                         \compute_QC_report.R
 
You can rename the `QC-dragdrop` folder to anything you like (try to avoid spaces - they usually cause trouble).
 
### Usage (short version)

  For the long version, see [user_manual.pdf].
  
  You can create a QC report by dragging a **txt-folder** (or any file within the txt-folder)
  onto the `createQC_dragNdrop.bat` file which resides on your (network) drive.


  [user_manual.pdf]: user_manual.pdf
  [http://www.r-project.org]: http://www.r-project.org