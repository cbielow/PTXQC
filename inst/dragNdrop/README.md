This file is for admins/support-staff who want to configure and install the PTXQC package
and provide an automated drag'n'drop solution to their MaxQuant users (incl. themselves) for creating QC reports.

**Note:** the drag'n'drop currently only works for Windows, not MacOS or Linux. 
          For the latter, please refer to the vignettes of this package (after installation), which explain how to invoke PTXQC from within R.

We recommend using a shared network drive as installation folder to which all users have read access. Alternatively you can copy the folder structure we are about to create to a folder on your local machine - however, if you want to use the QC reporting from multiple PC's just use a network folder.

When you are done, provide the users with the vignette `PTXQC-DragNDrop` (or write your own).
This will show them how to invoke the QC (it's very easy - see 'Usage' below).

### Installation
 - create clean R installation by downloading R from [http://www.r-project.org] and installing it. Make sure to install the x64 version for 64bit support.
   If you already have R installed, you might skip this step.
 - Now, we install PTXQC and its dependency packages.
   Start the R GUI and execute the following commands to install the `PTXQC` package (dependencies should be installed automatically):

   If you use an existing R installation which should not contain PTXQC, you can use the following snipped to create an extra library. It also avoids a copy of all your custom libraries unrelated to PTXQC to end up in the final installation.
   If you installed R anew or want PTXQC installed for your daily R work you should skip this section:
   
        ##
        ## the following section of code is only required to create a self-contained 
        ## library folder which hosts all required libraries for PTXQC for ease of copying
        ##
        lp = .libPaths() ## backup the names of the orignal lib paths
        tmp_dir = paste0(tempdir(), "\\tmp-PTXQC-Rpackages")
        dir.create(tmp_dir)
        .libPaths(tmp_dir)
        cat(paste0("\nNew temporary library folder: '", .libPaths()[1], "' created. Installing PTXQC ..."))

   Now, install PTXQC:
   
        ##
        ## the actual installation of packages
        ##
        if (!require(devtools, quietly = TRUE)) install.packages("devtools")
        library("devtools")             ## this might give a warning like 'WARNING: Rtools is required ...'. Ignore it.
        source("http://bioconductor.org/biocLite.R")
        biocLite("Biobase")
        install_github("cbielow/PTXQC", build_vignettes = TRUE) 
        help(package="PTXQC")           ## all done; check out the documentation

        cat(paste0("\nCopy the packages contained in '", .libPaths()[1], "' to your library folder (see installation manual)."))


 - The last command will tell you where R installed PTXQC. Open this folder (let's call it **_libR_**) in your file explorer. It can either be your default R library 
   folder which comes with R (e.g. `C:\Program Files\R\R-3.1.0\library`), or it could look like this `C:/Users/cbielow/AppData/Local/Temp/RtmporaVBo/tmp-PTXQC-Rpackages`. 
 - copy the folder **_libR_**`\PTXQC\inst\QC-dragdrop` to your target location (let's call it **_QCdir_**), usually some network or local drive
   (e.g. **_QCdir_** is `Z:\my-proteomics`)
 - copy the whole R installation directory (e.g. `c:\program files\R\R-3.1.0`) into the **_QCdir_**`\QC-dragdrop\_internal` sub-folder
 - rename your newly copied R installation directory to `R-3.1.0` (or edit the **_QCdir_**`\QC-dragdrop\createQC_dragNdrop.bat` to match your R version)
 
Now, you should have the following structure

    QCdir\QC-dragdrop\
                     \createQC_dragNdrop.bat
                     \createQC_dragNdrop_withYAML.bat
                     \_internal\
                               \R-3.1.0\...  (version may differ)
                               \compute_QC_report.R
 
You can rename the `QC-dragdrop` folder to anything you like (try to avoid spaces - they usually cause trouble).
 
### Usage (short version)

  For the long version, see the package vignette `PTXQC-DragNDrop`.
  
  You can create a QC report by dragging a **txt-folder** (or any file within the txt-folder)
  onto the `createQC_dragNdrop.bat` file which resides on your (network) drive.

  [http://www.r-project.org]: http://www.r-project.org