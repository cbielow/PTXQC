This file is for admins/support-staff who want to configure and install the PTXQC package
and provide an automated drag'n'drop solution to their MaxQuant users (incl. themselves) for creating QC reports.

All you need is a shared network drive to which all users have read access. Alternatively you can copy the folder structure we
are about to create to each user machine individually (but maintenance is going to be more work).

When you are done, provide the users with the [user_manual.pdf] from this folder (or write your own).
This will show them how to invoke the QC (it's very easy - see 'Usage' below).

### Installation
 - create clean R installation
 - install the PTXQC package (dependencies should be installed automatically)
 - if R installed the new packages to a custom location (e.g. your user home)
   , move them to the ./library within your new R installation folder
 - copy 'QC-dragdrop' onto a network drive (e.g. Z:\proteomics\QC-dragdrop)
 - copy the whole R installation directory (e.g. c:\program files\R-3.1.0) into the _internal(!) sub-folder (within 'QC-dragdrop')
 
Now, you should have the following structure

    <network-drive>\QC-dragdrop\
                               \createQC_dragNdrop.bat
                               \_internal\
                                         \R-3.1.0  (version may differ)
                                         \compute_QC_report.R
 
You can rename 'QC-dragdrop' to anything you like (try to avoid spaces - they usually cause trouble).
 
### Usage
  Now, you should be able to create a QC report by dragging a txt-folder (or any file within this folder)
  onto the createQC_dragNdrop.bat file which resides on your (network) drive.


  [user_manual.pdf]: user_manual.pdf