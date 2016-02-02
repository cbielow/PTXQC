This file is for admins/support-staff who want to configure and install the PTXQC package
and provide an automated drag'n'drop solution to their MaxQuant users (incl. themselves) for creating QC reports.
Please contact an admin or proficient R user, if you are new to installing R packages.

**Note:** the drag'n'drop currently only works for Windows, not MacOSX or Linux. 
          For the latter, please refer to the [vignettes] [Ref_Vign] of this package (after installation), which explain how to invoke PTXQC from within R.

We recommend using a shared network drive as installation folder to which all users have read access. Alternatively you can copy the folder structure we are about to create to a folder on your local machine. However, if you want to use the QC reporting from multiple PC's just use a network folder.

When you are done, provide the users with the [PTXQC-DragNDrop] [Ref_VignDrag] vignette (or write your own).
This will show them how to invoke the QC (it's very easy - see 'Usage' below).

### Installation
  
In order to make PTXQC available as a standalone application on a network drive, we will need to make a copy of the whole R installation
directory as one of the last steps. If you use an existing R installation which already contains a large package library, you might want to find a fresh PC without R
and install R from scratch, such that you do not end up copying useless packages which are not relevant for PTXQC.
If you are a knowledgeable R user already, you can probably find other ways around it (e.g. a temporary library, or selective copying of packages, but
this is out of scope of this manual).

 - if you want to build the vignettes which come with this package, you need to install [pandoc](https://github.com/jgm/pandoc/releases).
   Alternatively, you can skip vignette building and read the [vignettes] [Ref_Vign] online from the PTXQC GitHub page.
   **If you install pandoc later while your R session is already open, you need to close and re-open R in order to make R aware of pandoc!**

 - create a clean R installation by downloading R from [http://www.r-project.org] and installing it. Make sure to install the x64 version for 64bit support.
   If you already have R installed, you may skip this step.
 - Now, we install PTXQC and its dependency packages.
   Start the R GUI (64bit!) and execute the following commands to install the `PTXQC` package (dependencies are installed automatically).
   
   Run **each line** separately, i.e. do not copy and paste the whole block.
   If an error should occur, this allows to track it down more easily. See [FAQ - Installation] [Ref_VignFAQ]
   how to resolve them.
   If R is asking if it should create a user-specific library, say yes (this just means that R cannot write to it's global library, but that's fine).
   
        ##
        ## the actual installation of packages
        ##
        if (!require(devtools, quietly = TRUE)) install.packages("devtools")
        library("devtools")             ## if you see a warning like 'WARNING: Rtools is required ...': ignore it.
        install_github("cbielow/PTXQC", build_vignettes = TRUE)   ## use build_vignettes = FALSE if you did not install pandoc!

        cat(paste0("\nPTXQC was installed to '", .libPaths()[1], "'.\n\n"))

> **Here's some intuition about the next steps**   
> To make PTXQC Drag'n'drop a standalone application which can be run from any user and even from a network drive,
> we need a copy of the whole R installation. Within the PTXQC package (which part of the R installation) we have prepared a 
> batch script which invokes R and calls PTXQC. So first, we need to extract this script (i.e. copy it to a top-level folder)
> and then place the whole R installation beneath it. This includes all your R-libraries, which might reside in a different place.

Almost done, we just need to copy some folders.
        
 1. The last command will tell you where R installed PTXQC, i.e. `PTXQC was installed to '<libR>'.`
    **Open** this **`<libR>`** folder in your file explorer. It will either be your default R library 
    folder which comes with R (e.g. `C:\Program Files\R\R-3.1.0\library`), or a user-specific folder like 
    this `C:/Users/cbielow/Documents/R/win-library/3.1`. Keep this directory in mind. We will need it later.
 2. **Copy** the folder `<libR>\PTXQC\dragNdrop\QC-dragdrop` to a custom target location of your choice (let's call it **`<QCdir>`**) where you want PTXQC to reside.
    Usually, that's some network or local drive (e.g. **`<QCdir>`** = `Z:\my-proteomics`)
 3. **Copy** the whole R installation directory (e.g. `c:\program files\R\R-3.1.0`) into the `<QCdir>\QC-dragdrop\_internal` sub-folder,
    such that you end up with `<QCdir>\QC-dragdrop\_internal\R-3.1.0\` (your R version number might differ)
 4. If your newly copied R installation directory is _not_ named `R-3.1.0`, **rename** it to `R-3.1.0` 
    (or edit the `<QCdir>\QC-dragdrop\createQC_dragNdrop.bat` to match your R version)
 5. If the `<libR>` folder was not your default R library (see step 1), but a user library, then **copy** the packages it contains (PTXQC among others),
    into `<QCdir>\QC-dragdrop\_internal\R-3.1.0\library\`.
 
Now, you should have the following structure

    QCdir\QC-dragdrop\
                     \createQC_dragNdrop.bat
                     \createQC_dragNdrop_withYAML.bat
                     \_internal\
                               \R-3.1.0\...  (version may differ)
                               \compute_QC_report.R
 
You can rename the `QC-dragdrop` folder to anything you like (try to avoid spaces - they usually cause trouble).

To make PTXQC accessible more easily, you can create a shortcut on your Windows Desktop pointing to `createQC_dragNdrop.bat`. Use this shortcut as if it were the original batch file (see usage below).


### Updating the PTXQC Drag'n'drop installation

Start the R GUI within in PTXQC standalone installation at `<QCdir>\QC-dragdrop\_internal\R-3.1.0\bin\x64\rgui.exe` and
run the block of code from the installation (see above) again. Then, copy the packages from the `<libR>` folder to the main R library (see the last bullet point above).
 
 
### Usage (short version)

  For the long version, see the [PTXQC-DragNDrop] [Ref_VignDrag] vignette of this package.
  
  You can create a QC report by dragging a **txt-folder** (or any file within the txt-folder)
  onto the `createQC_dragNdrop.bat` file which resides on your (network) drive.


  [Ref_Vign]: https://github.com/cbielow/PTXQC/tree/master/vignettes
  [Ref_VignDrag]: https://github.com/cbielow/PTXQC/blob/master/vignettes/PTXQC-DragNDrop.Rmd
  [Ref_VignFAQ]: https://github.com/cbielow/PTXQC/blob/master/vignettes/PTXQC-FAQ.Rmd