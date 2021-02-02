
#'
#' Get the information of each CV term from an obo file.
#' 
#' @param cv_obo_file "xxx.obo"
#' @return A list containing cv term information
#' 
#' @export
#' 
parseOBO = function(cv_obo_file){
  ontology = ontologyIndex::get_ontology(cv_obo_file)
  obo = scan(file = cv_obo_file, what = "character")
  return(obo)
}


