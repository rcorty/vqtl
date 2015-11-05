#'  @title Get units of a scanonevar object
#'
#'  @author Robert Corty \email{rcorty@@gmail.com}
#'
#'
#'  @param object the scanonevar object whose units are interrogated
#'
#'  @return units of scanonevar object
#'
#'  @details none
units <- function(object) {

  return(attr(object, 'units'))
}
