qcMetric_EVD_modTable =  setRefClass(
  "qcMetric_EVD_modTable",
  contains = "qcMetric",
  methods = list(initialize=function() {  callSuper(    
    helpTextTemplate = 
      "Compute an occurence table of modifications (e.g. Oxidation (M)) for all peptides, including the unmodified.

The plot will show percentages, i.e. is normalized by the total number of peptide sequences (where different charge state counts as a separate peptide) per Raw file.

The sum of frequencies may exceed 100% per Raw file, since a peptide can have multiple modifications.
E.g. given three peptides in a single Raw file                 <br>
1. _M(Oxidation (M))LVLDEADEM(Oxidation (M))LNK_               <br>
2. _(Acetyl (Protein N-term))M(Oxidation (M))YGLLLENLSEYIK_    <br>
3. DPFIANGER                                                   <br>

, the following frequencies arise:

* 33% of 'Acetyl (Protein N-term)' <br>
* 33% of 'Oxidation (M)'           <br>
* 33% of '2 Oxidation (M)'         <br>
* 33% of 'Unmodified'              <br>

Thus, 33% of sequences are unmodified, implying 66% are modified at least once. 
If a modification, e.g. Oxidation(M), occurs multiple times in a single peptide it's listed as a separate modification (here '2 Oxidation (M)').

Heatmap score [EVD: Pep ModTable]: Deviation of (unmodified peptides fraction) when compared to a representative Raw file ('qualMedianDist' function).
",
    workerFcn = function(.self, df_evd)
    {
      ## completeness check
      if (!checkInput(c("fc.raw.file", "modifications"), df_evd)) return()

      name_unmod = "Unmodified"
      name_unmod_inverse = "Modified (total)"
      tbl = modsToTableByRaw(df_evd, name_unmod, name_unmod_inverse)
      lpl = byXflex(tbl, tbl$fc.raw.file, FUN = plot_peptideMods)
      #for (pl in lpl) print(pl)
      
      ## QC measure: deviation from representative sample (in terms of unmodified peptide count)
      qc_mods = tbl[tbl$modification_names %in% c(name_unmod, name_unmod_inverse),]
      qc_mods[, .self$qcName] = qualMedianDist(qc_mods$Freq / 100)
      
      return(list(plots = lpl, qcScores = qc_mods[, c("fc.raw.file", .self$qcName)]))
    }, 
    qcCat = "prep", 
    qcName = "EVD:~Peptide~VarMod", 
    orderNr = 0103
  )
    return(.self)
  })
)




