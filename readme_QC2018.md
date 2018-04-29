**Nutzung des Packages mit [mzTab](https://github.com/HUPO-PSI/mzTab) Files**

Zunächst muss das PTXQC Package geladen werden.
Da die von uns im Zuge des Praktikums bearbeitete Version nicht über den Cran-Server zu beziehen ist, muss das Github-Repository geladen werden.
Es wird der Branch [QC2018](https://github.com/cbielow/PTXQC/tree/QC2018) benötigt.

Um, falls nötig, das [githubinstall](https://cran.r-project.org/web/packages/githubinstall) Package und die Funktionen des QC2018 Branches des PTXQC Repositories zu laden, kann der folgende Befehl in die R Console kopiert werden.

    if (!require(githubinstall, quietly = TRUE)) install.packages("githubinstall")
    library("githubinstall")
    githubinstall("PTXQC", ref= "QC2018")
    
Folgender Befehl erstellt einen QC-Report anhand einer .mzTab:

    PTXQC::createReport("PATH/TO/FILE.mztab")
    
Es kann vorkommen, dass die Erstellung des Reports durch die [browser()](http://stat.ethz.ch/R-manual/R-devel/library/base/html/browser.html) Funktion unterbrochen wird. In diesem Fall kann der Vorgang durch Enter oder Continue fortgesetzt werden.
Die browser() Funktion ist hilfreich zum Debugging. Mit ihr können zB. Variablenbelegungen während der Laufzeit überprüft werden.



